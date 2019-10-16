import os
os.chdir(r'C:\Sabine\fasta_sabine\SVs_Fasta_perseus\20190329_extractFASTAMQ')

import fastaFunctions
import pandas as pd

import re

accessions, sequences = fastaFunctions.read_fasta_file_fullID('uniprot_proteome_rattus_norvergicus.fasta')

UniProtIDs = []

uniprotPattern = re.compile(r'.+\|(\S+)\|.+')

for acc in accessions:
    if uniprotPattern.match(acc):
        match = uniprotPattern.match(acc)
        thisID = match.group(1)
        UniProtIDs.append(thisID)
    

ID_srt = pd.read_excel(r'C:\Sabine\fasta_sabine\SVs_Fasta_perseus\ID_for_fasta_20190408.xlsx', sheetname='ID_for_fasta_sortiert').iloc[0:400,:]

ID_all = pd.read_excel(r'C:\Sabine\fasta_sabine\SVs_Fasta_perseus\ID_for_fasta_20190408.xlsx', sheetname='ID_for_fasta_20190408').iloc[0:400,:]


# proteinGroups = proteinGroups[proteinGroups['Potential contaminant'] != '+']
# proteinGroups = proteinGroups[proteinGroups['Reverse'] != '+']

# proteinGroups.sort_values(by='iBAQ', ascending=False, inplace=True)


# write fasta for first entrys in majority protein IDs
intenseIDs = ID_srt['Majority protein IDs'].values.tolist()

mostIntenseProteins = []

for ID in intenseIDs:
    temp = ID.split(';')
    # mostIntenseProteins.extend(temp)
    mostIntenseProteins.append(temp[0])

with open('ID_srt_filtered.fasta', 'w') as out:
    for protein in mostIntenseProteins:
        print('Looking for {}'.format(protein))
        found = False
        for idx, thisID in enumerate(UniProtIDs):
                if protein == thisID:
                    print('Found {}'.format(thisID))
                    out.write('>' + accessions[idx] + '\n')
                    out.write(sequences[idx] + '\n')
                    print('Wrote {}'.format(accessions[idx]))
                    found = True
        if found == False:
            print('Couldnt find {} in DB'.format(protein))


# write fasta for first entrys in majority protein IDs
intenseIDs = ID_all['Majority protein IDs'].values.tolist()

mostIntenseProteins = []

for ID in intenseIDs:
    temp = ID.split(';')
    mostIntenseProteins.extend(temp)
    # mostIntenseProteins.append(temp[0])

with open('ID_all_filtered.fasta', 'w') as out:
    for protein in mostIntenseProteins:
        found = False
        for idx, thisID in enumerate(UniProtIDs):
                if protein == thisID:
                    out.write('>' + accessions[idx] + '\n')
                    out.write(sequences[idx] + '\n')
                    found = True
        if found == False:
            print('Couldnt find {} in DB'.format(protein))