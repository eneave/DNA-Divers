################################################
## Dowloading UK fish reference library from  ##
################################################

## CODE AND REFLIB FROM https://github.com/genner-lab/meta-fish-lib

### START A FRESH R SESSION ###

# load packages (install if required)
library("tidyverse")
library("ape")

# load REMOTE references and cleaning scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

# subset reference library table by metabarcode fragment (primer set) from the following options:
print(tibble(metabarcodes=c("coi.lerayxt","coi.ward","12s.miya","12s.riaz","12s.valentini","12s.taberlet","16s.berry","cytb.minamoto")))
# change 'metabarcode' argument as appropriate:
reflib.sub1 <- subset_references(df=reflib.cleaned, metabarcode="12s.taberlet")
reflib.sub2 <- subset_references(df=reflib.cleaned, metabarcode="12s.miya")

# [OPTIONAL] taxonomically dereplicate and filter on sequence length
# 'proplen=0.5' removes sequences shorter than 50% of median sequence length
# 'proplen=0' retains all sequences
reflib.sub1 <- derep_filter(df=reflib.sub1, derep=TRUE, proplen=0.5)
reflib.sub2 <- derep_filter(df=reflib.sub2, derep=TRUE, proplen=0.5)


# write out reference library in FASTA and CSV format to current working directory
# currently supported fasta formats are: [1] sintax, [2] dada2 (taxonomy), [3] dada2 (species), and [4] plain dbid (GenBank or BOLD database identifiers)
write_references_fasta(df=reflib.sub1)
write_references_fasta(df=reflib.sub2)







