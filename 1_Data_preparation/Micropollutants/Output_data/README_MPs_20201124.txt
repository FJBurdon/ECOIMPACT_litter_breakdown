ECOIMPACT Toxcicity data - MPs converted to Toxic Units
README - Dr. Francis J. Burdon

1. A .csv file including max values for sum TUs for HM, PB, and Other TUs has been subsetted from the data used in the second macroinvertebrate paper (Burdon et al. 2017). This is based on based on tox data from all invertebrates, using the geometric mean for critical values. Values for U1 are used as dummy variables for U2 - this is for data exploration purposes and not final analyses (use only D and U1 locations):

-       Output: MP_Ecoimpact_invert_TU_final.csv
-	Input: Tab_Delta_allMacro_v20181116.csv

2. A .csv file has similarly been generated from the data file containing median values for sum TUs for HM, PB, and Other (also Insecticides, Herbicides, and Fungicides). The max values for the three broad groups (HM, PB, and Other) are provided above. This data includes tox data from all invertebrates, and Dapnia magna (DM) using the geometric mean and logarithmic mean. These data correspond to the data provided by Nicole Munz. Values for U1 are used as dummy variables for U2 - this is for data exploration purposes and not final analyses (use only D and U1 locations):

-	Output: MP_Ecoimpact_sumTU_msPAF_final.csv
-	Input: Spear_median_sumTU_msPAF.csv 

The above input data file can be found the switchdrive (switchdrive\EcoImpact_Macro\Output Data\Overview_sumTU_msPAF\spear_median_sumTU_msPAF.csv).

3. A .csv file has similarly been generated from the raw data file used by Nicole Munz to generate the above data using critical toxicity values for Dapnia magna. This data file includes mean, median, and max values for sum TUs: total, HM, Insecticides, Herbicides, Fungicides, Pharmaceuticals, and Other. This data includes tox data from Dapnia magna (DM) using the geometric mean. These data correspond to the data provided by Nicole Munz. Values for U1 are used as dummy variables for U2 - this is for data exploration purposes and not final analyses (use only D and U1 locations):

-	Output: 3_DATA_MP_TUs_metasubstances_final.csv
-	Input: 3_DATA_MPs_TU_metasubstance_VER02_161114.csv 

The above output data (and output from 1 and 2) are generated in the R script "Rscript_MPs_data_preparation_201124.R" 
