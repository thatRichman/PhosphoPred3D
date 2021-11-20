import urllib.parse
import urllib.request
import pandas as pd
import itertools

from pandas.core.groupby import grouper

# reclosedev stackoverflow
def grouper_it(n, iterable):
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)

psp_df = pd.read_csv("/home/spencer/PhosphoPred3D/data/Phosphorylation_site_dataset_psp.csv")
up_accs = psp_df.ACC_ID.unique().tolist()

url = "https://www.uniprot.org/uploadlists/"
with open("../data/PDB_IDS.tab", 'w+') as out:
    out.write("UP_ACC\tPDB_ID\n")
    for chunk in grouper_it(100, up_accs):
        params = {
            'from': 'ACC+ID',
            'to': 'PDB_ID',
            'format': 'tab',
            'query': " ".join(chunk)
            }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            resp = f.read().decode('utf-8').split("\n")[1:]
            out_dat = [x + "\n" for x in resp if x is not None]
        out.writelines(out_dat)