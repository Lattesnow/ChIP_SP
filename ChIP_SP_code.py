import pandas as pd

# input the hic matrix file from Mustache
hic_df = pd.read_excel('*.tsv')

# input the ChIP-bed peak file from MACS2 call peak 
chipseq_df = pd.read_excel('chipseq_bed_data.bed')

# Create an empty list to store the matching rows
matching_rows = []

# Loop through each row of the ChIP-seq bed data
for idx, row in chipseq_df.iterrows():
    # Filter the hic matrix data to only include rows that match the current ChIP-seq row's chromosome
    filtered_hic_df = hic_df[hic_df['BIN1 CHR'] == row['chr']]
    # Filter the filtered_hic_df to only include rows that overlap with the ChIP-seq region
    filtered_hic_df = filtered_hic_df[((filtered_hic_df['BIN1 START'] <= row['end']) & (filtered_hic_df['BIN1 END'] >= row['end'])) |
                                      ((filtered_hic_df['BIN1 START'] <= row['start']) & (filtered_hic_df['BIN1 END'] >= row['start'])) |
                                      ((filtered_hic_df['BIN1 START'] <= row['start']) & (filtered_hic_df['BIN1 END'] >= row['end']))]
    # Check if any rows matched the criteria and append them to the list of matching rows
    if not filtered_hic_df.empty:
        for _, hic_row in filtered_hic_df.iterrows():
            matching_rows.append([hic_row['BIN1 CHR'], hic_row['BIN1 START'], hic_row['BIN1 END'], hic_row['BIN2 CHR'], hic_row['BIN2 START'], hic_row['BIN2 END'], row['start'], row['end'], row['score'], row['strand']])

# Convert the list of matching rows to a dataframe and write it to an Excel file
output_df = pd.DataFrame(matching_rows, columns=['BIN1 CHR', 'BIN1 START', 'BIN1 END', 'BIN2 CHR', 'BIN2 START', 'BIN2 END', 'start', 'end', 'score', 'strand'])
output_df.to_excel('output.xls', index=False)
