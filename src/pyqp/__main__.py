"""Quantitative Proteomic Service

Usage:
    pyqp api
    pyqp cli <proteomicTSV> <proteomeXML> [--field=<quantity_column>] [--adress=<apiAdress>] [--port=<apiPort>] [--verbose] [--topScore=<pathway_number>]

Options:
  -h --help     Show this screen.
  --field=<quantity column>  csv column header featuring signal 
  --purb=purb  aa
  --intg=intg  bbb
  --alpha=alpha  ccc
  --ncore=ncore  ddd
  --sizelim=sizelim eee
  --prot=<proteomeXML>  ggg
  --adress=<apiAdress>  aaa
  --port=<apiPort>  aaa
  --verbose  iiii
  --topScore=<pathway_number>  aaaa

"""

# TEST W/ mycoplasma proteome
# The test this
#python -m pyqp cli previous/wt2_subset.tsv unigo/src/unigo/data/uniprot-proteome_UP000000625.xml.gz

from docopt import docopt

#from pyT2GA import analysis
from unigo import Unigo as createUniGOTree
from unigo import uloads as createGOTreeFromAPI

from .utils import proteomicWrapper
from pyproteinsExt.uniprot import EntrySet as createUniprotCollection

from requests import get
from .api import app
import time



arguments = docopt(__doc__)

#print(arguments)
abnd_field = arguments['--field'] if arguments['--field'] else "Corrected Abundance ratio (1,526968203)"
nTop = int(arguments['--topScore']) if arguments['--topScore'] else 20

if arguments['cli']:
    quantProteomic = proteomicWrapper(            csv_file = arguments['<proteomicTSV>'], abnd_label = abnd_field)
    uColl          = createUniprotCollection(collectionXML = arguments['<proteomeXML>']  )


    missingProt = []
    for x in quantProteomic.uniprot:
        if not uColl.has(x):
            print(f"{x} not found in proteome")
            missingProt.append(x)
        
    for x in missingProt:
        quantProteomic.remove(x)

  
    taxid = uColl.taxids[0]
    apiAdress = arguments['--adress'] if arguments['--adress'] else "127.0.0.1"
    apiPort   = arguments['--port']   if arguments['--port']   else "5000"
    url = f"http://{apiAdress}:{apiPort}/unigo/{taxid}"
    print(f"Fetching universal annotation tree from {url}")
    
    expUniprotID = [ _ for _ in quantProteomic.uniprot ]
    resp = get(url)
    if resp.status_code == 404:
        print(f"{url} returned 404, provided proteome XML {taxid} may not be registred")
    else:
        unigoTree = createGOTreeFromAPI(resp.text, expUniprotID)
        x,y = unigoTree.dimensions
        print("Unigo Object successfully buildt w/ following dimensions:")
        print(f"\txpTree => nodes:{x[0]} children_links:{x[1]}, total_protein_occurences:{x[2]}, protein_set:{x[3]}")  
        print(f"\t universeTree => nodes:{y[0]} children_links:{y[1]}, total_protein_occurences:{y[2]}, protein_set:{y[3]}")  
        
        nDelta=int(0.1 * len(quantProteomic))
        print(f"{len(quantProteomic)} proteins available in quantitative records, taking first {nDelta} as of quantity modified")
        print("Computing ORA")
        deltaUniprotID = expUniprotID[:nDelta]
        rankingsORA = unigoTree.computeORA(deltaUniprotID, verbose = arguments['--verbose'])
        print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")    
    

if arguments['api']:
    app.run(port=1234)

"""    
    unigoTree = createUniGOTree( backgroundUniColl = uColl,
                                proteinList       = [ x for x in quantProteomic.uniprot ],
                                fetchLatest       = False)


    start = time.perf_counter()
    # Taking 10% w/ highest qtty value
    rankingsORA = unigoTree.computeORA(
                        [ _ for _ in quantProteomic[nTop].uniprot ]
                        , verbose = False)
    stop = time.perf_counter()
    print(f"Test Top - {5}\n{rankingsORA[5]}")
    print(f"Execution time is {stop-start} sc")
"""

# Unnecssary
def typeGuardTaxID(proteomicData, uColl):
    taxids = {}
    for uID in proteomicData.uniprot:
        uObj = uColl.get(uID)
        if not uObj.taxid in taxids:
            taxids[uObj.taxid] = 0
        taxids[uObj.taxid] += 1
    return sorted( [ (k,v) for k,v in taxids.items() ], key=lambda x:x[1] )

#r = pyt2ga.analysis(proteoRes, GOpwRes, STRINGRes, mapperRes, intg=False,
#                   abnd_label = "Corrected Abundance ratio (1,526968203)", ncore=3)