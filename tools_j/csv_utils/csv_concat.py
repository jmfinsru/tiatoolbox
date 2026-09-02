import pandas as pd

"""
Combining two CSV files that contain the cell counts for the 40x samples from inference.
"""
#CSV file paths
csv1 = '/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/Func116_ST_HE_40x_BF_01/counts/nuclei_counts_from_0_to_349.csv'
csv2 = '/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/Func116_ST_HE_40x_BF_01/counts/nuclei_counts_from_351_to_699.csv'
csv3 = '/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/Func116_ST_HE_40x_BF_01/counts/nuclei_counts_from_701_to_1087.csv'
df1 = pd.read_csv(csv1) #Load CSV files
df2 = pd.read_csv(csv2)
df3 = pd.read_csv(csv3)
combined = pd.concat([df1, df2, df3], ignore_index=True) #Combine CSV files

#Save as CSV file
combined.to_csv('/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/Func116_ST_HE_40x_BF_01/counts/nuclei_counts_combined.csv', index=False)