from pysam import VariantFile
import pandas as pd
import csv
import numpy as np

vcf_filename = "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
panel_filename = "igsr-1000 genomes 30x on grch38.tsv.tsv"

# get a dictionary with sample IDs and population codes
with open(panel_filename) as panel_file:
    labels = {}  # {sample_id: population_code}
    panel_file.readline()
    for line in panel_file:
        line = line.strip().split('\t')
        labels[line[0]] = line[3]
    
population_data = pd.read_csv(panel_filename, sep = '\t')

ancestries = population_data['Superpopulation name'].unique()
print(ancestries)

asian_superpops = ['South Asia (SGDP),South Asian Ancestry',
                 'East Asia (SGDP),East Asian Ancestry', 
                 'East Asian Ancestry', 'South Asian Ancestry',
                 'European Ancestry,West Eurasia (SGDP)']
columns = ['Sample name', 'Population code', 'Superpopulation name']
df = population_data[columns]
asian_df = df[df['Superpopulation name'].isin(asian_superpops)]
                    
print(asian_df)

asian_samples = list(asian_df['Sample name'])
asian_samples_df = pd.DataFrame(asian_samples)
# now we have a df with asian samples and their corresponding population codes 

print(asian_samples_df)
asian_samples_df.to_csv('asian_samples.csv', index = False, header = False)

vcf = VariantFile(vcf_filename)  # num_rows = 1066557

# now, time to filter alleles by sample ID
sample_ids = list(vcf.header.samples)

##dff = pd.read_csv('asian_samples.csv', header = None)


filtered_vcf_filename = "filtered-vcf.vcf"
filtered_vcf = VariantFile(filtered_vcf_filename)
filtered_sample_ids = list(filtered_vcf.header.samples)

#extracting genotypes from the filtered data
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
    if counter % 10666 == 0:
        print(counter)
        print(f'{round(100 * counter / 1066557)}%')

#print(variant_ids)
genotypes = np.array(genotypes)
print(genotypes.shape)

matrix = np.count_nonzero(genotypes, axis=2)

matrix = matrix.T
print(matrix.shape) 

# filter the other sample IDs and population codes out of the labels dictionary

filtered_labels = {k:v for k,v in labels.items() if k in filtered_sample_ids}
print(len(filtered_labels))

print(len(sample_ids))


final_df = pd.DataFrame(matrix, columns = variant_ids, index = sample_ids)
final_df['Population code'] = final_df.index.map(filtered_labels)
final_df.to_csv("asian_population_matrix.csv")

# dict_you_want = { your_key: old_dict[your_key] for your_key in your_keys }
# foodict = {k: v for k, v in mydict.items() if k.startswith('foo')}


