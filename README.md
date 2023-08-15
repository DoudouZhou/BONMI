# BONMI
Codes and data for the paper ``Multi-source learning via completion of block-wise overlapping noisy matrices'' (https://arxiv.org/abs/2105.10360).


## Data
In the paper, we used four different sources:

- 20 million notes at Stanford (Finlayson et al., 2014). The link to the preprocessed PMI matrix: https://www.dropbox.com/s/hcd50pf5z7pufau/ppmi.stan.concepts_perBin_30d.RData
-  10 million notes of 62K patients at Partners Healthcare System (PHS): the data is not shareable.
-  Health records from MIMIC-III, a freely accessible critical care database (Johnson et al., 2016): the raw data can be obtained from https://physionet.org/content/mimiciii/1.4/ but we are not allowed to distribute it. 
-  A Chinese PMI matrix from multiple sources including medical textbooks and Wikipedia: https://www.dropbox.com/s/9cxn8u8ahspt4xm/PPMI_Other8_window_10threshold_0weighted_0.rds



# Reference

- Finlayson, S. G., LePendu, P., and Shah, N. H. (2014). Building the graph of medicine from millions of clinical narratives. Scientific data, 1:140032.

- Johnson, A. E., Pollard, T. J., Shen, L., Li-Wei, H. L., Feng, M., Ghassemi, M., Moody, B., Szolovits, P., Celi, L. A., and Mark, R. G. (2016). Mimic-iii, a freely accessible critical care database. Scientific data, 3(1):1–9.
