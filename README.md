# CRLM_analysis
Aim of this computational pipeline is to combine and explore pathological annotation and data from clinical reports on colorectal cancer liver metastases (e.g. regression after chemotherapy, growth pattern distribution, survival in subgroups).

An example dataset is provided in the input directory. It contains randomly generated clinical and pathological report data as well as random annotation ndpa files. These data are completely artificial and do not reflect any patient data. Invasion front annotations were done with NDP.view 2 and are exported as .ndpa files into the folder "./input/annotations". Further details are provided in the manuscript. 

The pipeline is divided in several parts, whose function can be summarized as following: 

1_ reads the ndpa files and converts each .ndpa file into one dataframe which is saved as .csv files in "./output/parsed_annotations"

2_ combines all annotations into one dataframe which is saved as .csv file in "./output/combined_annotations.csv"

3_ tidies the annotations into growth pattern annotations per probe/slide/tumor and tumor regression per probe/slide/tumor and saves them in "./output/"

4_ imports two clinical dataframes from .csv files, reorganizes the data, calculates the age of the patients, plots an age histogram, merges the clinical with pathological report data, calculates tumor regression, calculates tumor diameter sums, imports growth pattern annotation, reorganizes and categorizes growth pattern annotations

5_ explores growth pattern abundance in the CRLM dataset, (i) per probe, (ii) per tumor, (iii) per slide in all patients and grouped for chemonaive versus neoadjuvant chemotherapy

6_ explores the degree of tumor regression in the CRLM dataset

7_ explores the overall and liver progression-free survival in the CRLM dataset, depending on various factors, and performs  uni- as well as multivariate analyses

9_ builds a clinical data table with (i) all patients combined (ii) patients grouped for chemonaive and neoadjuvant chemotherapy 

10_ correlates the degree of desmoplastic growth with tumor regression per probe with (i) all patients combined (ii) patients grouped for chemonaive and neoadjuvant chemotherapy 

To further show the functionality of the pipeline 4_, 5_, 6_, 7_, 9_, 10_ are provided as html rendered R markdown files.
