#Parallel Reaction Monitoring Script
#Written by Sameer Varma and Stanley M. Stevens Jr.
#Last updated: May 25, 2023
#Written with Pandas and tested with Python3
#Numpy is used only for formatting floating point output
#Beta 0.5 fixes:
# 1) MS filetype supported: DIANN and TIMSDIANN
#Beta 0.4 fixes:
# 1) Account for multiple entries in Genes (to activate, set split_genes = 'YES')
# 2) Account for incomplete digestion
# 3) Account for variability in sequence charge states
# 4) Account for modified amino acids (UniMod)
# 5) Output control -- top N intensity sequences + drop sequences with certain amino acids.

import pandas as pd
import numpy as np
import re
import os
import datetime

debug = 'OFF' # or ON


def split(word):
    return [char for char in word]

def assign_mass(mysequence,mycharge):
    #monoisotopic masses taken from https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
    masses = {'A':71.03711, 'L':113.08406, 'V':99.06841, 'I':113.08406, 'Y':163.06333, 'P':97.05276, 'F':147.06841, 'W':186.07931, 'D':115.02694, 'E':129.04259, 'R':156.10111, 'K':128.09496, 'N':114.04293, 'Q':128.05858, 'C':103.00919, 'S':87.03203, 'T':101.04768, 'G':57.02146, 'M':131.04049, 'H':137.05891}
    mass_hydrogen=1.00727647
    mass_water=18.01057
    if debug == 'ON': print("assigning masses to unmodified sequences")
    if debug == 'ON': print(mysequence,mycharge)
    split_sequence = split(mysequence)
    if debug == 'ON': print(split_sequence)
    mymass=0
    for v in split_sequence:
        mymass+=masses.get(v);
    mymass+=mycharge*mass_hydrogen
    mymass+=mass_water
    if debug == 'ON': print(mymass)
    return(mymass)

def adjust_UniMod(mysequence):
    #monoisotopic masses taken from https://unimod.org
    UniMod_masses = {1:42.010565, 4:57.021464, 7:0.984016, 35:15.994915}
    if debug == 'ON': print("assigning mass correction to account for amino acid modifications")
    mymass=0
    mymods = [int(x) for x in re.findall(r'\d+', mysequence)]
    if debug == 'ON': print(mymods)
    for x in mymods:
        mymass+=UniMod_masses.get(x);
    if debug == 'ON': print(mymass)
    return(mymass)

def number_KR(mysequence):
    count=0
    split_sequence = split(mysequence)
    for v in split_sequence:
        if v=='R': count+=1
        if v=='K': count+=1
    return(count)


#INITIALIZE variables

#'Genes' in df_MS can have multiple entries separated by semicolon
#To split 'Genes' entries and create identical copies for each 'Genes' entries, set split_genes to YES
split_genes = 'NO'

#Isolation Width definitions
Isolation_width_cutoff = 700
Isolation_width_low = 2
Isolation_width_high= 3

#Tolerances
Intensity_charge_state_tolerance = 0.2 # 20%
Intensity_digestion_tolerance = 0.1 # 10%
Min_fragment_length = 0
RT_range = 240
IM_tolerance = 4

#MS file types supported
supported_MS_filetypes = ["DIANN", "TIMSDIANN"]
#default
filetype_MS = 'DIANN'

#Output control
filename_out = 'output.csv' #default file name
topN = -1 #Set to a positive integer to print only the top N highest intensity sequences. To print all sequences, set topN = -1
drop_sequence_aaM = 'NO' #Set to YES to drop sequences consisting of Methionine
mydecimals = 5 #for formatting output


print("----------------------------------------------")
print("****  Parallel Reaction Monitoring (PRM)  ****")
print("****      Version: Beta 0.5 (05/2023)     ****")
print("****          Stevens group @ USF         ****")
print("----------------------------------------------")
currentTime = datetime.datetime.now()
if currentTime.hour < 12: print('Good morning',os.getlogin(),"! \n")
elif 12 <= currentTime.hour < 18: print('Good afternoon',os.getlogin(),"! \n")
else: print('Good evening',os.getlogin(),"! \n")

print("Defaults:")
print("Debug mode:", debug)
print("Isolation width cutoff:",Isolation_width_cutoff)
print("Isolation width low:",Isolation_width_low)
print("Isolation width high:",Isolation_width_high)
print("Charge-state relative intensity tolerance:",Intensity_charge_state_tolerance*100,"%")
print("Undigested sequence relative intensity tolerance:",Intensity_digestion_tolerance*100,"%")
print("Smallest fragment length:",Min_fragment_length)
print("Split genes:", split_genes)

#Read IO file names
print('-'*80)
print("Input prompts:")
filename_IPA = input("Enter name of IPA file (txt format): ")
filename_MS = input("Enter name of MS file (tsv format): ")
while True:
    filetype_MS = input("Enter MS file type (TIMSDIANN/DIANN):")
    if filetype_MS.upper() in supported_MS_filetypes:
       break
    else:
        print('Unsupported filetype, please enter again:')
filename_out = input("Enter name of output file (csv format): ")
RT_range = float(input("Enter RT range (seconds): [Suggested:240] "))
IM_tolerance = float(input("Enter IM tolerance (percent): [Suggested:4] "))
Intensity_charge_state_tolerance = float(input("Enter charge-state relative intensity tolerance (fraction): [Suggested:0.2] "))
Intensity_digestion_tolerance = float(input("Enter undigested sequence relative intensity tolerance (fraction): [Suggested:0.1] "))
topN = int(input("Enter a positive integer to print only the top N highest intensity sequences. (Enter -1 for all sequences):"))
drop_sequence_aaM = input("Drop seqeunces contaning Methionine (YES/NO): ")
print('-'*80+'\n')

#TESTING DONE USING THE FOLLOWING
#filename_IPA = 'Hela_AutophagyIPA.txt'
#filename_MS  = 'DIANNout_HelaDIA.tsv'
#filetype_MS  = 'DIANN'
#filename_MS  = 'TIMSDIANNout_HelaDIA.tsv'
#filetype_MS  = 'TIMSDIANN'
#filename_out = 'test.csv'

#Output columns
if debug == 'ON':
    column_names = ['Mass_over_Charge','Precursor.Charge','Isolation_width','RT_median','RT_range','Start_IM','End_IM','CE','External_ID','Description','Intensity_mean','Digestion', 'Modified.Sequence']
    #Write header to output file
    with open(filename_out, 'w+') as file:
        file.write('Mass [m/z],Charge,Isolation Width [m/z],RT [s],RT Range [s],Start IM [1/K0],End IM [1/K0],CE [eV],External ID,Description,MaxLFQ_mean,Digestion\n')
else:
    column_names = ['Mass_over_Charge','Precursor.Charge','Isolation_width','RT_median','RT_range','Start_IM','End_IM','CE','External_ID','Description']
    #Write header to output file
    with open(filename_out, 'w+') as file:
        file.write('Mass [m/z],Charge,Isolation Width [m/z],RT [s],RT Range [s],Start IM [1/K0],End IM [1/K0],CE [eV],External ID,Description\n')


print("\nAnalyzing Input Files")

#Assign dataFrames for input files
df_IPA = pd.read_csv(filename_IPA,delimiter="\t",skiprows=2)
if debug == 'ON': print(df_IPA)
print("Genes in IPA file:",len(df_IPA))
df_MS = pd.read_csv(filename_MS,delimiter="\t")
if debug == 'ON': print(df_MS)
print("Rows in MS file:",len(df_MS))
#'Genes' in df_MS can have multiple entries separated by semicolon
#Split 'Genes' entries by creating identical copies for each entry
if split_genes == 'YES':
    df_MS['Genes'] = df_MS['Genes'].str.split(';')
    df_MS = df_MS.explode('Genes').reset_index(drop=True)
    cols = list(df_MS.columns)
    cols.append(cols.pop(cols.index('File.Name')))
    df_MS = df_MS[cols]
    print("Adding rows in MS file to accomodate multiple entries in 'Genes' column:",len(df_MS))

#Convert strings to lowercase for comparison
df_IPA['Symbol'] = df_IPA['Symbol'].str.lower()
df_MS['Genes'] = df_MS['Genes'].str.lower()

#Iterating over List of genes in IPA
print("\nGrouping...")
for g_ipa in df_IPA['Symbol']:
    #Find symbol entries in genes
    df_MS_select = df_MS.loc[df_MS['Genes']==g_ipa]
    #if group not empty
    print('-'*50)
    print("Searching for Gene",g_ipa,": found",len(df_MS_select),"entries")
    print('-'*50)
    if(len(df_MS_select)==0):
        print("---------------WARNING: ZERO ENTRIES FOR GENE",g_ipa,"---------------") 
    else:
        #Calculate averages, etc for each group of sequences and add to dataframe
        if filetype_MS.upper()== 'DIANN':
            df_MS_select_grouped = df_MS_select.groupby(['Modified.Sequence','Precursor.Charge']).agg({'RT':['median'],'IM':['median'], 'Precursor.Normalised':['mean'],'Stripped.Sequence':['first']})
        else:
            df_MS_select_grouped = df_MS_select.groupby(['Modified.Sequence','Precursor.Charge']).agg({'RT':['median'],'Exp.1/K0':['median'], 'Precursor.Normalised':['mean'],'Stripped.Sequence':['first']})
        df_MS_select_grouped.columns = ['RT_median','IM_median','Intensity_mean','Stripped.Sequence']
        df_MS_select_grouped = df_MS_select_grouped.reset_index()
        df_MS_select_grouped['Gene'] = g_ipa
        df_MS_select_grouped['RT_median'] = 60*df_MS_select_grouped['RT_median']
        df_MS_select_grouped['RT_median'] = np.round(df_MS_select_grouped['RT_median'], decimals=mydecimals)
        df_MS_select_grouped['RT_range']  = RT_range
        df_MS_select_grouped['Start_IM'] = df_MS_select_grouped['IM_median']*(1 - IM_tolerance/100)
        df_MS_select_grouped['Start_IM']=np.round(df_MS_select_grouped['Start_IM'], decimals=mydecimals)
        df_MS_select_grouped['End_IM'] = df_MS_select_grouped['IM_median']*(1 + IM_tolerance/100)
        df_MS_select_grouped['End_IM']=np.round(df_MS_select_grouped['End_IM'], decimals=mydecimals)
        df_MS_select_grouped['Intensity_mean'] = np.round(df_MS_select_grouped['Intensity_mean'], decimals=mydecimals)
        
        #Assign masses and isolation widths and create description
        for i in df_MS_select_grouped.index:
            df_MS_select_grouped.at[i,'Description'] = df_MS_select_grouped.at[i,'Gene'] + '_' + df_MS_select_grouped.at[i,'Modified.Sequence']
            #Calculate mass from raw sequence
            mass=assign_mass(df_MS_select_grouped.at[i,'Stripped.Sequence'],df_MS_select_grouped.at[i,'Precursor.Charge'])
            #Adjust for chemical modifications of amino acids
            mass+=adjust_UniMod(df_MS_select_grouped.at[i,'Modified.Sequence'])
            #calculate mass/charge
            df_MS_select_grouped.at[i,'Mass_over_Charge'] = mass/df_MS_select_grouped.at[i,'Precursor.Charge']
            #assign isolation width
            if(df_MS_select_grouped.at[i,'Mass_over_Charge']<Isolation_width_cutoff):
                df_MS_select_grouped.at[i,'Isolation_width']=Isolation_width_low
            else:
                df_MS_select_grouped.at[i,'Isolation_width']=Isolation_width_high
        df_MS_select_grouped['Mass_over_Charge']=np.round(df_MS_select_grouped['Mass_over_Charge'], decimals=mydecimals)

        if debug == 'ON': print("Grouped according to Precursor.Charge and Modified.Sequence")
        if debug == 'ON': print(df_MS_select_grouped)
      
        #Add empty coloumns
        df_MS_select_grouped['CE']=""
        df_MS_select_grouped['External_ID']=""

        #Examine intensities of identitcal sequences with different charge states
        print("Examining sequences with different charges")
        indices_dropped=[]
        for i in df_MS_select_grouped.index:
            for j in df_MS_select_grouped.index:
                if (df_MS_select_grouped.at[j,'Stripped.Sequence'] == df_MS_select_grouped.at[i,'Stripped.Sequence']) and (j != i):
                    Intensity_sum = df_MS_select_grouped.at[i,'Intensity_mean'] + df_MS_select_grouped.at[j,'Intensity_mean']
                    Intensity_j_rel = df_MS_select_grouped.at[j,'Intensity_mean']/Intensity_sum
                    if Intensity_j_rel < Intensity_charge_state_tolerance:
                        print(g_ipa,": Removing sequence", df_MS_select_grouped.at[j,'Stripped.Sequence'], "with charge", df_MS_select_grouped.at[j,'Precursor.Charge'], "as its Intensity=",df_MS_select_grouped.at[j,'Intensity_mean'], "is less than", Intensity_charge_state_tolerance, "of total Intensity of", Intensity_sum)
                        indices_dropped.append(j)
        indices_dropped = [*set(indices_dropped)]
        df_MS_select_grouped = df_MS_select_grouped.drop(labels=indices_dropped)
        print("Total sequences dropped = ", len(indices_dropped))
        if debug == 'ON': print("Indices of sequences dropped:", indices_dropped)

        #Look for incomplete digestion
        print("Examining incomplete digestion")
        indices_dropped=[]
        for i in df_MS_select_grouped.index:
            my_number_KR = number_KR(df_MS_select_grouped.at[i,'Stripped.Sequence'])
            if my_number_KR>1:
                df_MS_select_grouped.at[i,'Digestion'] = 'Incomplete'
                for j in df_MS_select_grouped.index:
                    if (i != j):
                        if (df_MS_select_grouped.at[j,'Stripped.Sequence'] in df_MS_select_grouped.at[i,'Stripped.Sequence']) and (df_MS_select_grouped.at[i,'Stripped.Sequence'] != df_MS_select_grouped.at[j,'Stripped.Sequence']) and (len(df_MS_select_grouped.at[i,'Stripped.Sequence'])-len(df_MS_select_grouped.at[j,'Stripped.Sequence'])>Min_fragment_length):
                            if (df_MS_select_grouped.at[i,'Intensity_mean'] < Intensity_digestion_tolerance * df_MS_select_grouped.at[j,'Intensity_mean']):
                                indices_dropped.append(i)
                                print("Dropping undigested fragment",df_MS_select_grouped.at[i,'Stripped.Sequence'], "as its Intensity=", df_MS_select_grouped.at[i,'Intensity_mean'], "is", Intensity_digestion_tolerance, "of the Intensity of",  df_MS_select_grouped.at[j,'Intensity_mean'], "of its digested fragment", df_MS_select_grouped.at[j,'Stripped.Sequence'])
                            else:
                                indices_dropped.append(i)
                                indices_dropped.append(j)
                                print("Dropping undigested fragment",df_MS_select_grouped.at[i,'Stripped.Sequence'], "and its digested fragment", df_MS_select_grouped.at[j,'Stripped.Sequence'], "as their Intensities of", df_MS_select_grouped.at[i,'Intensity_mean'], "and", df_MS_select_grouped.at[j,'Intensity_mean'], "are within", Intensity_digestion_tolerance, "of each other." )
                
            else:
                df_MS_select_grouped.at[i,'Digestion'] = 'Complete'
        indices_dropped = [*set(indices_dropped)]
        df_MS_select_grouped = df_MS_select_grouped.drop(labels=indices_dropped)
        print("Total sequences dropped = ", len(indices_dropped))
        if debug == 'ON': print("Indices of sequences dropped:", indices_dropped)
            
                #Split sequence
                #fragments=re.split(r'R|K',df_MS_select_grouped.at[i, 'Stripped.Sequence'])
                #fragments=[x for x in fragments if x != '']
                #print(fragments, len(fragments))
            
        #drop sequences containing Methionine
        if drop_sequence_aaM == 'YES':
            indices_dropped=[]
            for i in df_MS_select_grouped.index:
                if ('M' in df_MS_select_grouped.at[i,'Stripped.Sequence']): indices_dropped.append(i)
            df_MS_select_grouped = df_MS_select_grouped.drop(labels=indices_dropped)
            print("Sequences containing methionine dropped = ", len(indices_dropped))
            if debug == 'ON': print("Indices of sequences dropped:", indices_dropped)
            
        #Sort by Intensity
        df_sorted = df_MS_select_grouped.sort_values(by=['Intensity_mean'],ascending=False)
        if debug == 'ON': print(df_sorted)

        #Write group to OUTPUT file
        if topN == -1:
            df_sorted.to_csv(filename_out, mode="a", index=False, header=False, columns=column_names)
        else:
            df_get_rows = df_sorted.head(topN)
            if debug == 'ON': print(df_get_rows)
            df_get_rows.to_csv(filename_out, mode="a", index=False, header=False, columns=column_names)

print("\nNORMAL TERMINATION: Finished writing output file", filename_out, "\n")
