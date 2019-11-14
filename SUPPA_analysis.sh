#!/usr/bin/bash

SUPPA analysis of AS
============================================

### generate events from annotation
python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py generateEvents -i AtRTD2_QUASI.gtf -e SE SS RI FL MX -o AtRTD2_all_events.ioe -f ioe

## make tx tpm dataframe from kallisto counts 
python3 ~/.local/lib/python3.4/site-packages/SUPPA/multipleFieldSelection.py -i in/kallisto/*/abundance.tsv -k 1 -f 5 -o out/iso_all_sample.tpm.txt

## calculate psiPerEvent
python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py psiPerEvent -i in/AtRTD2_all_events.ioe -e out/iso_all_sample.tpm.txt -o out/psy_allsamples


### differential splicing analysis

  A=OE_ASCO
  B=Col

# split isoform quantification table  
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/iso_all_sample.tpm.txt X35S.ASCO_1_S7_noribo_clean.kallisto,X35S.ASCO_2_S8_noribo_clean.kallisto,X35S.ASCO_3_S9_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.tpm out/$B.tpm  -i

# split psi quantification table
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/psy_allsamples.psi X35S.ASCO_1_S7_noribo_clean.kallisto,X35S.ASCO_2_S8_noribo_clean.kallisto,X35S.ASCO_3_S9_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.psi out/$B.psi  -e
  
# calculate differential splicing
  mkdir out/${A}_vs_${B}
  
  python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py diffSplice -m empirical -gc -i AtRTD2_all_events.ioe \
  -p out/${A}.psi out/${B}.psi -e out/${A}.tpm out/${B}.tpm -o out/${A}_vs_${B}/comp
  
########################################
 
  A=RNAi_ASCO
  B=Col
# split isoform quantification table 
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/iso_all_sample.tpm.txt RNAi.ASCO_1_S4_noribo_clean.kallisto,RNAi.ASCO_2_S5_noribo_clean.kallisto,RNAi.ASCO_3_S6_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.tpm out/$B.tpm  -i
 
# split psi quantification table
  A=RNAi_ASCO
  B=Col
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/psy_allsamples.psi RNAi.ASCO_1_S4_noribo_clean.kallisto,RNAi.ASCO_2_S5_noribo_clean.kallisto,RNAi.ASCO_3_S6_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.psi out/$B.psi  -e
  
# calculate differential splicing
  mkdir out/${A}_vs_${B}
  
  python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py diffSplice -m empirical -gc -i AtRTD2_all_events.ioe \
  -p out/${A}.psi out/${B}.psi -e out/${A}.tpm out/${B}.tpm -o out/${A}_vs_${B}/comp
