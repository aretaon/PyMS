# -*- coding: utf-8 -*-
"""
Skript downloads a set of fasta-files using the EBI dbfetch tool
starting from a list of Uniprot identifiers.

See: https://www.biostars.org/p/90271/
"""

import re
import urllib.request
import math

pattern = re.compile(r'(P{1}\d{5})')

accessions = []

txt_path = r'Z:\02_experiments\04_mona_copII_microsomes\04c_cop_interaction_bs3\20180405_jb04c002_generate_fasta\Partnerproteine.txt'

out_path = r'Z:\02_experiments\04_mona_copII_microsomes\04c_cop_interaction_bs3\20180405_jb04c002_generate_fasta\Partnerproteine'

with open(txt_path, 'r') as f:
    for line in f.readlines():
        if pattern.match(line):
            acc = pattern.match(line).group(1)
            accessions.append(acc)

# remove duplicates
accessions = list(set(accessions))

for i in range(0,math.ceil(len(accessions)/200)):
    
    # slice a set of < 200 entries from the DB as max 200 are allowed from EBI dbfetch
    curr_accessions = accessions[i*200: i*200 + 200]
    
    # Download the fasta using the EBI dbfetch tool
    outpath = out_path + '_{}.fasta'.format(i+1)

    pre = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id='
    
    identifiers = '+'.join(curr_accessions)
    
    post = '&format=fasta&style=raw&Retrieve=Retrieve'
    
    # urllib works as wget request do retrieve and download the fasta file from
    # the EBI server
    urllib.request.urlretrieve(pre+identifiers+post, outpath)

outcontent = []
for i in range(0,math.ceil(len(accessions)/200)):
    fasta_file = outpath = out_path + '_{}.fasta'.format(i+1)
    with open(fasta_file, 'r') as f:
        for line in f.readlines():
            outcontent.append(line)

with open(out_path + '.fasta', 'w') as o:
    o.writelines(outcontent)