import pandas as pd
import numpy as np
from pysam import VariantFile

panel_filename = "igsr-1000 genomes 30x on grch38.tsv.tsv"

# obtaining a df with asian samples and 
# their corresponding population codes and
# superpopulation names 
population_data = pd.read_csv(panel_filename, sep = '\t')
asian_superpops = ['South Asia (SGDP),South Asian Ancestry',
                 'East Asia (SGDP),East Asian Ancestry', 
                 'East Asian Ancestry', 'South Asian Ancestry',
                 'European Ancestry,West Eurasia (SGDP)']
columns = ['Sample name', 'Population code', 'Superpopulation name']
df = population_data[columns]
asian_df = df[df['Superpopulation name'].isin(asian_superpops)]
                    
print(asian_df) 
 
# saving sample IDs to a csv file
# used to filter alleles from asian samples in the vcf file 
asian_samples = list(asian_df['Sample name'])
asian_samples_df = pd.DataFrame(asian_samples)
asian_samples_df.to_csv('asian_samples.csv', index = False, header = False)


# extracting genotypes from the filtered data
filtered_vcf_filename = "filtered-vcf.vcf" 
filtered_vcf = VariantFile(filtered_vcf_filename)

genotypes = []
variant_ids = []
sample_ids = []
counter = 0
for record in filtered_vcf:
    counter += 1
    if counter % 100 == 0: 
        alleles = [record.samples[x].allele_indices for x in record.samples]
        sample_ids = [sample for sample in record.samples]
        genotypes.append(alleles)
        variant_ids.append(record.id)
    if counter % 10666 == 0:  # num_rows = 1066557
        print(counter)
        print(f'{round(100 * counter / 1066557)}%')


genotypes = np.array(genotypes)
print(genotypes.shape)

matrix = np.count_nonzero(genotypes, axis=2)

matrix = matrix.T
print(matrix.shape) 

# get a dictionary with sample IDs and population codes
with open(panel_filename) as panel_file:
    labels = {}  # {sample_id: population_code}
    panel_file.readline()
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[3]

# filtering non-asian sample IDs and 
# population codes out of the labels dictionary
filtered_sample_ids = list(filtered_vcf.header.samples)
filtered_labels = {k:v for k,v in labels.items() if k in filtered_sample_ids}

# exporting allele, sample ID, and population code data
final_df = pd.DataFrame(matrix, columns = variant_ids, index = sample_ids)
final_df['Population code'] = final_df.index.map(filtered_labels)
final_df.to_csv("asian_population_matrix.csv")