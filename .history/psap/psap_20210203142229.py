"""Main module."""
import os
import sys
import ntpath
import datetime
import pandas as pd
from Bio import SeqIO
from scipy import signal
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm.auto import tqdm
from multiprocessing import Pool
import time
from pathlib import Path

uniprot_ids = '''Q99700
O60885
Q14781
P45973
P38432
Q7Z5Q1
O00571
Q9NQI0
O43781
Q04637
Q15056
P15502
Q01844
P22087
Q06787
P35637
Q13283
Q9UN86
P10071
P62993
Q13151
P09651
Q32P51
P22626
P51991
O14979
P31943
P55795
P31942
O43561
P10636
P43243
Q15648
P19338
P06748
Q15233
P52948
Q01860
P11940
P29590
Q8WXF1
Q96PK6
P98179
P23246
Q16637
P00441
Q07889
P23497
Q07955
Q01130
O95793
O75683
Q92804
Q13148
Q15554
P31483
Q01085
Q9UHD9
P46937
Q15059
P10644
P54727
Q13501
Q9NPI6
O00444
Q53HL2
Q9ULW0
Q9Y6A5
Q96LT7
P06748
Q16236
O43823
Q07157
Q9UDY2
O95049
Q9Y6M1
Q6NT89
Q8NC56
P04150
Q9BYJ9
Q9Y5A9
Q7Z739
P10997
Q9GZV5
P51608
O14979
P43351
Q08379
P08621
Q15424
Q16630
P18615
P48443
'''

RESIDUES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Kyte & Doolittle {kd} index of hydrophobicity
HP = {'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5,
      'Q':-3.5, 'E':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5,
      'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6,
      'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2, 'U': 0.0}

class MakeMatrix:
    def __init__(self, dbfasta):  
        self.df = pd.DataFrame()
        executables = [
             'self.fasta2df(dbfasta)',
#             'self.amino_acid_analysis()',
#             'self.idr_iupred()',
#             'self.hydrophobic()',
#             'self.add_iupred_features()',
#             'self.add_hydrophobic_features()',
#             'self.add_biochemical_combinations()',
#             'self.add_lowcomplexity_features()' ,
#             #'self.add_plaac()'
        ]        
        for e in executables:
            start = time.time()     
            print(e)
            exec(e)
            end = time.time()        
            print(str(round(end - start, 2))+'s '+e)

    def fasta2df(self, dbfasta):
        rows = list()
        with open(dbfasta) as f:
            for record in SeqIO.parse(dbfasta, 'fasta'):
                seqdict = dict()
                seq = str(record.seq)
                id = record.description.split('|')                
                if id[0] == 'sp':
                    uniprot_id = id[1]
                    name = id[2].split(' ')[0]
                    rows.append([name, uniprot_id, seq])
                elif id[0] == 'tr':
                    uniprot_id = id[1]
                    name = id[2].split(' ')[0]
                    rows.append([name, uniprot_id, seq])
                else:
                    uniprot_id = id[0]
                    name = id[2].split(' ')[0]
                    rows.append([name, uniprot_id, seq])                    
        self.df = pd.DataFrame(rows, columns=['protein_name', 'uniprot_id', 'sequence'])
        
def export_matrix(name, fasta_path, out_path, operating_system='Windows', annotate_uniprot_ids=True):
    # Change pathing
    """ Generates and saves a file which contains features of a protein sequence.
    Parameters:
        name: Name of the file.
        fasta_path: Path of the fasta file which needs to be featured.
        operating_system: String which indicates which operating system is used only 'Windows' available.
    """
    data = MakeMatrix(fasta_path)   
    now = datetime.datetime.now()
    date = (str(now.day) + '-' + str(now.month)  + '-' +  str(now.year))
    #if operating_system == 'Windows':
    pkl = Path(out_path,name+'_llps_f2f_'+date+'.pkl')
    data.df.to_pickle(pkl)
    print(pkl)
    if ()
    
    
    return data.df

def annotate(df, uniprot_ids, identifier_name):
    uniprot_ids = [s.strip() for s in uniprot_ids.splitlines()]
    df[identifier_name] = 0
    for prot_id in uniprot_ids:
        df.loc[df['uniprot_id'] == prot_id, identifier_name] = 1
        if (len(df.loc[df['uniprot_id'] == prot_id])) == 0:
            print(prot_id+' is not found.')
    return df