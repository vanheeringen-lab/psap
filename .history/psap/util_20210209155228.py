from psap import MakeMatrix
import datetime
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
    date = (str(now.day)+'-'+str(now.month)+'-'+str(now.year))
    #if operating_system == 'Windows':
    pkl = Path(out_path, name+'_'+identifier_name+'_'+date+'.pkl')
    data.df.to_pickle(pkl)
    print(pkl)
    print("Adding labels to df")
    df_ann = annotate(data.df, uniprot_ids, identifier_name)
    pkl_ann = Path(out_path,name+'_'+identifier_name+'_'+date+'_ann_'+'.pkl')
    df_ann.to_pickle(pkl_ann)
    print(pkl_ann)