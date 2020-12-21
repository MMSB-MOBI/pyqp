"""Quantitative Proteomic Service

Usage:
    pyqp web
    pyqp cli <proteomicTSV> <proteomeXML> [--field=<quantity column>]

Options:
  -h --help     Show this screen.
  --field=<quantity column>  csv column header featuring signal 
  --purb=purb  aa
  --intg=intg  bbb
  --alpha=alpha  ccc
  --ncore=ncore  ddd
  --sizelim=sizelim eee
  
"""

# TEST W/ mycoplasma proteome
# The test this
#python -m pyqp cli previous/wt2_subset.tsv unigo/src/unigo/data/uniprot-proteome_UP000000625.xml.gz

from docopt import docopt

#from pyT2GA import analysis
from unigo import Unigo as createUniGOTree
from .utils import proteomicWrapper
from pyproteinsExt.uniprot import EntrySet as createUniprotCollection
arguments = docopt(__doc__)

print(arguments)
abnd_field = arguments['--field'] if arguments['--field'] else "Corrected Abundance ratio (1,526968203)"

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

    nTop=int(0.1 * len(quantProteomic))
    print(f"{len(quantProteomic)} proteins available in quatitative records, taking first {nTop}")

    unigoTree = createUniGOTree( backgroundUniColl = uColl,
                                proteinList       = [ x for x in quantProteomic.uniprot ],
                                fetchLatest       = False)


    # Taking 10% w/ highest qtty value
    rankingsORA = unigoTree.computeORA(
                        [ _ for _ in quantProteomic[nTop].uniprot ]
                        , verbose = False)
    
    print(f"Test Top - {5}\n{rankingsORA[5]}")




#r = pyt2ga.analysis(proteoRes, GOpwRes, STRINGRes, mapperRes, intg=False,
#                   abnd_label = "Corrected Abundance ratio (1,526968203)", ncore=3)