cd /home/agatawa/archetypes/omics

for I in *.prot.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_prot
# 29
for I in *.myriad.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_myriad
## 12
for I in *.OLINK.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_olink
## 28
for I in *.targ.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_targ_met
# 22 
for I in *.untarg.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_untarg_met
# 29
for I in *.transcript.txt; do awk '{if ($5>0 && $4<0.01){print $0}}' $I | sort -k 4 -g | head; done | awk '{print $1}' | sort -u  > top_genes
## 27

### find total unique proteins differing between subgroups: total=490
for I in *.prot.txt; do awk '{if ($4<0.05){print $1}}' $I ; done | sort -u  > all_signif_proteins
for I in *.myriad.txt; do awk '{if ($4<0.05){print $1}}' $I ; done | sort -u  >> all_signif_proteins
for I in *.OLINK.txt; do awk '{if ($4<0.05){print $1}}' $I ; done | sort -u  >> all_signif_proteins

### find total unique metabolites differing between subgroups: total=277
for I in *.targ.txt; do awk '{if ($4<0.05){print $1}}' $I ; done | sort -u  > all_signif_metabs
for I in *.untarg.txt; do awk '{if ($4<0.05){print $1}}' $I ; done | sort -u  >> all_signif_metabs

