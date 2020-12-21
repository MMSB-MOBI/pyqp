import pandas as pd

from unigo import Unigo as createGOTree



class proteomicWrapper():
    def __init__(self, csv_file=None, df=None, abnd_label=None):
        if csv_file:
            self.pdFrame = pd.read_csv(csv_file, sep='\t')
        else:
            self.pdFrame = df.reset_index(drop=True)

        if abnd_label:
            if not abnd_label in self.pdFrame.columns:
                raise ValueError(f"column {abnd_label} not found in tabular data")

            self.pdFrame = self.pdFrame.sort_values(by=[abnd_label], ascending=False).dropna(subset=[abnd_label])
            self.pdFrame = self.pdFrame.reset_index(drop=True)

    def __getitem__(self, _slice):
        return proteomicWrapper(df=self.pdFrame[_slice])
        
    def head(self):
        pass
    @property 
    def uniprot(self):
        for index, row  in self.pdFrame.iterrows():
            yield row['Accession']        
        
    def remove(self, uniprotID):
        df = self.pdFrame
        df.drop(df.loc[df['Accession']==uniprotID].index, axis=0, inplace=True)
    def __len__(self):
        return len(self.pdFrame.index)

def create_STRING(uniColl, proteomicData, proxy=None, file=None): # Do we implement a file ??
    """

        http_proxy  = "http://ftprox.ibcp.fr:3128"
        https_proxy = "http://ftprox.ibcp.fr:3128"
        ftp_proxy   = "http://ftprox.ibcp.fr:3128"
        proxy = { 
                  "http"  : http_proxy, 
                  "https" : https_proxy, 
                  "ftp"   : ftp_proxy
                }
    """
    
    tagSTRING  = "protein.physical.links.v11.0"
    # Ensure STRING_taxid is constant accross protein collection
    # Generate uniprotID<->STRING_ID mapper
    uniprotSTRING_mapper = [("String_id", "Uniprot_id")]
    _uniprotID = None
    
    for uniprotID in proteomicData.uniprot:
        e = uniColl.get(uniprotID)
        if e.STRING_ID:
            taxID, protID =  e.STRING_ID.split('.')
            if not _uniprotID is None:
                if _taxID != taxID:
                    print(f"Warning uniprot entries features different STRING taxon identifiers [{_uniprotID}]{_taxID} {_protID} != [{uniprotID}]{taxID} {protID} <== Droping this one")
                    continue
            _uniprotID, _taxID, _protID = uniprotID, taxID, protID
            uniprotSTRING_mapper.append( (protID.replace(f"{taxID}.",""), uniprotID) )
        else:
            print(f"Warning No STRING ID available for {uniprotID}")
    taxid = taxID

    # Fetching taxon specific string and format it
    headREST_url = "https://stringdb-static.org/download"
    tail_REST_url = f"{tagSTRING}.txt.gz"

    url=f"{headREST_url}/{tagSTRING}/{taxid}.{tail_REST_url}"
    print(url)

    r = requests.get(url, proxies=proxy)
    file = gzip.open(io.BytesIO(r.content))

    data=[("protein1","protein2", "experimental")]
    file.readline()
    for l in file.readlines():
            _ = l.decode('UTF-8').rstrip().split()
            data.append(  (_[0].replace(f"{taxid}.",""), _[1].replace(f"{taxid}.",""), float(_[2])/1000 )  )
    
    return data, uniprotSTRING_mapper

# DONOT USE
def create_GO_pathways(uniprotCollection, proteomicData,
                       root    = "biological_process",
                       owlFile = None
                      ):
    GOtree = createGOTree(
                   backgroundUniColl = uniprotCollection,
                   proteinList       = [ x for x in proteomicData.uniprot ],
                   owlFile = owlFile
                   )
    GO_pws = {}

    ok=False
    for node in GOtree.walk():
        if ok:    
            GO_pws[node.ID] = {
                "Name"     : node.name, 
                "Proteins" : node.getMembers(nr=True)
            }
        if node.name == root:
            ok = True
    return GO_pws, GOtree