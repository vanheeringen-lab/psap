"""Main module."""
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm.auto import tqdm
import time
from scipy import signal

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
HP = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
      'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
      'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
      'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2, 'U': 0.0}


class MakeMatrix:
    def __init__(self, dbfasta):
        self.df = pd.DataFrame()
        self.add_features()

    def add_features(self):
        executables = [
             'self.fasta2df(dbfasta)',
             'self.amino_acid_analysis()',
             'self.add_biochemical_combinations()',
             'self.add_lowcomplexity_features()'
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

    def hydrophobic(self):
        for index, row in self.df.iterrows():
            hpilst = pd.Series(list(row['sequence'])).map(HP).tolist()
            self.df.loc[index, 'HydroPhobicIndex'] = HydroPhobicIndex(hpilst)
            

    def amino_acid_analysis(self):
        for res in RESIDUES:
            self.df['fraction_'+res] = self.df['sequence'].str.count(res) / self.df['sequence'].str.len()
        self.df['length'] = self.df['sequence'].str.len()
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
        #for index, row in self.df.iterrows():
            seq = row['sequence']
            seqanalysis = ProteinAnalysis(seq)
            acidist = seqanalysis.get_amino_acids_percent()
            self.df.loc[index, 'IEP'] = seqanalysis.isoelectric_point()
            if 'X' not in seq and 'B' not in seq:
                self.df.loc[index, 'molecular_weight'] = seqanalysis.molecular_weight()
            if 'U' not in seq and 'X' not in seq and 'B' not in seq:
                self.df.loc[index, 'gravy'] = seqanalysis.gravy()

    def add_iupred_features(self):
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
        #for index, row in self.df.iterrows():
            idr = row['iupred'].glob[0]
            self.df.loc[index, 'idr_percetage'] = sum(i > .5 for i in list(idr))
            self.df.loc[index, 'idr_50'] = sum(i > .5 for i in list(idr)) / len(str(row['sequence']))
            self.df.loc[index, 'idr_60'] = sum(i > .6 for i in list(idr)) / len(str(row['sequence']))
            self.df.loc[index, 'idr_70'] = sum(i > .7 for i in list(idr)) / len(str(row['sequence']))
            self.df.loc[index, 'idr_80'] = sum(i > .8 for i in list(idr)) / len(str(row['sequence']))
            self.df.loc[index, 'idr_90'] = sum(i > .9 for i in list(idr)) / len(str(row['sequence']))
       
    @staticmethod   
    def convolve_signal(sig, window=25):
        win = signal.hann(window)
        sig = signal.convolve(sig, win, mode='same') / sum(win)
        return sig

    
    def add_hydrophobic_features(self):
        hpi0, hpi1, hpi2, hpi3, hpi4, hpi5 = list(), list(), list(), list(), list(), list()
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
        #for index, row in self.df.iterrows():
            sw = convolve_signal(row['HydroPhobicIndex'].hpilist, window=30)
            hpi0.append(sum(i < -1.5 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-1.5_frac'] = hpi
            hpi1.append(sum(i < -2.0 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-2.0_frac'] = hpi
            hpi2.append(sum(i < -2.5 for i in sw) / len(sw))
            # self.df.loc[index, 'hpi_<-2.5_frac'] = hpi
            hpi3.append(sum(i < -1.5 for i in sw))
            # self.df.loc[index, 'hpi_<-1.5'] = hpi
            hpi4.append( sum(i < -2.0 for i in sw))
            # self.df.loc[index, 'hpi_<-2.0'] = hpi
            hpi5.append(sum(i < -2.5 for i in sw))
            # self.df.loc[index, 'hpi_<-2.5'] = hpi
        self.df['hpi_<-1.5_frac'] = hpi0
        self.df['hpi_<-2.0_frac'] = hpi1
        self.df['hpi_<-2.5_frac'] = hpi2
        self.df['hpi_<-1.5'] = hpi3
        self.df['hpi_<-2.0'] = hpi4
        self.df['hpi_<-2.5'] = hpi5

    def add_biochemical_combinations(self):
        df = self.df
        df = df.assign(Asx=df['fraction_D'] + df['fraction_N'])
        df = df.assign(Glx=df['fraction_E'] + df['fraction_Q'])
        df = df.assign(Xle=df['fraction_I'] + df['fraction_L'])
        df = df.assign(Pos_charge=df['fraction_K'] + df['fraction_R'] + df['fraction_H'])
        df = df.assign(Neg_charge=df['fraction_D'] + df['fraction_E'])
        df = df.assign(Aromatic=df['fraction_F'] + df['fraction_W'] + df['fraction_Y'] + df['fraction_H'])
        df = df.assign(Alipatic=df['fraction_V'] + df['fraction_I'] + df['fraction_L'] + df['fraction_M'])
        df = df.assign(Small=df['fraction_P'] + df['fraction_G'] + df['fraction_A'] + df['fraction_S'])
        df = df.assign(Hydrophilic=(df['fraction_S'] + df['fraction_T'] + df['fraction_H'] +
                                    df['fraction_N'] + df['fraction_Q'] + df['fraction_E'] +
                                    df['fraction_D'] + df['fraction_K'] + df['fraction_R']))
        df = df.assign(Hydrophobic= (df['fraction_V'] + df['fraction_I'] + df['fraction_L'] +
                                     df['fraction_F'] + df['fraction_W'] + df['fraction_Y'] +
                                     df['fraction_M']))
        # Added in version 2
        for dimer in ['GV', 'VG', 'VP', 'PG', 'FG', 'RG', 'GR', 'GG', 'YG', 'GS', 'SG', 'GA', 'GF', 'GD', 'DS']:
            self.df[dimer] = self.df['sequence'].str.count(dimer)
        df = df.assign(alpha_helix=df['fraction_V'] + df['fraction_I'] + df['fraction_Y'] + df['fraction_F']
                      + df['fraction_W'] + df['fraction_L'])
        df = df.assign(beta_turn=df['fraction_N'] + df['fraction_P'] + df['fraction_G'] + df['fraction_S'])
        df = df.assign(beta_sheet=df['fraction_E'] + df['fraction_M'] + df['fraction_A'] + df['fraction_L'])
        #Calculates the aromaticity value of a protein according to Lobry, 1994.
        # It is simply the relative frequency of Phe+Trp+Tyr.
        df = df.assign(aromaticity=df['fraction_F'] + df['fraction_W'] + df['fraction_Y'])
        self.df = df
        del df

    def add_lowcomplexityscore(self):
        lcs_window = 20
        lcs_cutoff = 7
        for index, row in self.df.iterrows():
            seq = str(row['sequence'])
            if len(seq) > lcs_window+1:
                sig = list()
                for i in range(len(seq)):
                    window = (seq[i: i+lcs_window])
                    if len(window) == lcs_window:
                        acid_comp = len(list(set(window)))
                        sig.append(acid_comp)
                score = sum([1 if i <= 7 else 0 for i in sig])
                self.df.loc[index, 'lcs_score'] = score
                self.df.loc[index, 'lcs_fraction'] = score / len(sig)

    def add_lowcomplexity_features(self):
        n_window = 20
        cutoff = 7
        n_halfwindow = int(n_window / 2)
        lcs_lowest_complexity = list()
        lcs_scores = list()
        lcs_fractions = list()
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0]):
        #for index, row in self.df.iterrows():
            # Determine low complexity scores
            seq = str(row['sequence'])
            lcs_acids = list()
            sig = list()
            # New
            lc_bool = [False] * len(seq)
            for i in range(len(seq)):
                if i < n_halfwindow:
                    peptide = seq[:n_window]
                elif i+n_halfwindow > int(len(seq)):
                    peptide = seq[-n_window:]
                else:
                    peptide = seq[i-n_halfwindow:i+n_halfwindow]
                complexity = (len(set(peptide)))
                if complexity <= 7:
                    for bool_index in (i-n_halfwindow, i+n_halfwindow):
                        try:
                            lc_bool[bool_index] = True
                        except IndexError:
                            pass
                    lcs_acids.append(seq[i])
                sig.append(complexity)
            # Adding low complexity scores to list
            low_complexity_list = pd.DataFrame({'bool':lc_bool, 'acid':list(seq)}, index=None)
            lcs_lowest_complexity.append(min(sig))
            lcs_scores.append(len(low_complexity_list.loc[low_complexity_list['bool'] == True]))
            lcs_fractions.append(len(low_complexity_list.loc[low_complexity_list['bool'] == True]) / len(seq))
            low_complexity_list = pd.DataFrame({'bool':lc_bool, 'acid':list(seq)}, index=None)
            if len(lcs_acids) >= n_window:
                for i in RESIDUES:
                    self.df.loc[index ,i+'_lcscore'] = (len(low_complexity_list.loc[
                        (low_complexity_list['bool'] == True) &
                        (low_complexity_list['acid'] == i)])
                    )
                    self.df.loc[index ,i+'_lcfraction'] = (len(low_complexity_list.loc[
                        (low_complexity_list['bool'] == True) &
                        (low_complexity_list['acid'] == i)]) / len(lcs_acids)
                    )
        self.df['lcs_fractions'] = lcs_fractions
        self.df['lcs_scores'] = lcs_scores
        self.df['lcs_lowest_complexity'] = lcs_lowest_complexity

    def add_plaac(self):
        plaac = CWD+'/data/plaac/plaac_swissprot140219.tsv'
        plaac = pd.read_csv(plaac, sep='\t')
        plaac[['database', 'accession', 'name']] = plaac['SEQid'].str.split('|',expand=True)
        plaac = plaac.drop(['SEQid', 'database', 'name', 'PAPAaa', 'STARTaa',
                            'ENDaa', 'COREaa', 'MW', 'MWstart', 'MWend', 'MWlen'], axis=1)
        self.df = pd.merge(self.df, plaac, left_on='uniprot_id', right_on='accession')
        self.df = self.df.drop('accession', axis=1)


def convolve_signal(sig, window=25):
    win = signal.hann(window)
    sig = signal.convolve(sig, win, mode='same') / sum(win)
    return sig


def annotate(df, uniprot_ids, identifier_name):
    uniprot_ids = [s.strip() for s in uniprot_ids.splitlines()]
    df[identifier_name] = 0
    for prot_id in uniprot_ids:
        df.loc[df['uniprot_id'] == prot_id, identifier_name] = 1
        if (len(df.loc[df['uniprot_id'] == prot_id])) == 0:
            print(prot_id+' is not found.')
    return df


def export_matrix(name, fasta_path, out_path, identifier_name=""):
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
    pkl = Path(out_path,name+'_'+identifier_name+'_'+date+'.pkl')
    data.df.to_pickle(pkl)
    print(pkl)
    print("Adding labels to df")
    df_ann = annotate(data.df, uniprot_ids, identifier_name)
    pkl_ann = Path(out_path,name+'_'+identifier_name+'_'+date+'_ann_'+'.pkl')
    df_ann.to_pickle(pkl_ann)
    print(pkl_ann)



