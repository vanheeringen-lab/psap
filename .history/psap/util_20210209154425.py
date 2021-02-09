from scipy import signal
from  psap import MakeMatrix

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