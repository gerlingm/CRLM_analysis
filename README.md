# CRLM_analysis
Combine and explore pathological annotation and clinical data on colorectal cancer liver metastases (e.g. regression after chemotherapy, growth pattern distribution)

Invasion front annotations are exported as .ndpa files into the folder "./input/annotations"

1_ reads the ndpa files and converts each .ndpa file into one dataframe which is saved as .csv files in "./output/parsed_annotations"

2_ combines all annotations into one dataframe which is saved as .csv file in "./output/combined_annotations.csv"

3_ tidies the annotations into growth pattern annotations per probe/slide/tumor and tumor regression per probe/slide/tumor and saves them in "./output/"

4_ imports two clinical dataframes from .csv files, reorganizes the data, calculates the age of the patients, plots an age histogram, merges the clinical with pathological report data
