/usr/bin/bash

SUPPA analysis gene of AS
============================================


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
