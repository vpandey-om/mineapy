import pandas as pd 
import os 

def getErythromycin_data():   
    # Define the path to the Excel file
    file_path = os.path.join('..',  'input', 'spectrum.00317-23-s0002.txt')

    # Read the Excel file into a DataFrame
    df = pd.read_csv(file_path,sep='\t', skiprows=1)
    
    # Specify the columns you want to extract
    columns_to_extract = ['Gene_id', 'ERY_group_fpkm', 'log2FoldChange(ERYvsETH)', 'pval(ERYvsETH)', 'padj(ERYvsETH)', 'significant(ERYvsETH)']

    # Extract the specified columns from the DataFrame
    df_selected = df[columns_to_extract]

    # Display the first few rows of the extracted DataFrame
    return df_selected 

def getLiverStage():
    # Define the path to the Excel file
    #file_path = os.path.join('..',  'input', 'male_female.txt')
    #file_path = os.path.join(os.path.dirname(__file__), 'input', 'malaria_life_cycle_bulk_gene_expression_data.txt')
    file_path = os.path.join('..',  'input', 'malaria_life_cycle_bulk_gene_expression_data.txt')
    # print(file_path)
    # Read the Excel file into a DataFrame
    df = pd.read_csv(file_path,sep='\t')

    ## gene length path
   
    #df_gene_length = pd.read_csv(os.path.join(os.path.dirname(__file__),  'input', 'GeneLength_Summary.txt'),sep='\t')
    df_gene_length = pd.read_csv(os.path.join('..',  'input', 'GeneLength_Summary.txt'),sep='\t')
    # Merge the DataFrames on the gene identifier columns
    merged_df = pd.merge(df, df_gene_length, left_on='Unnamed: 0', right_on='Gene ID', how='left')
    
    # Assuming `merged_df` is your DataFrame and has the required columns

    # Remove rows where 'Transcript Length' is NaN or 0 to avoid division errors
    filtered_df = merged_df.dropna(subset=['Transcript Length'])
    filtered_df = filtered_df[filtered_df['Transcript Length'] > 0]

    # Calculate total counts (sum of expression values) for each column
    total_counts_A = filtered_df['EEF_54h_A'].sum() / 1e6  # Scaling to millions
    total_counts_B = filtered_df['EEF_54h_B'].sum() / 1e6  # Scaling to millions

    # Calculate RPKM for each column
    filtered_df['RPKM_EEF_54h_A'] = (
        filtered_df['EEF_54h_A'] * 1e9 / (filtered_df['Transcript Length'] * total_counts_A)
    )

    filtered_df['RPKM_EEF_54h_B'] = (
        filtered_df['EEF_54h_B'] * 1e9 / (filtered_df['Transcript Length'] * total_counts_B)
    )

    # Calculate the mean of the RPKM values from 'EEF_54h_A' and 'EEF_54h_B'
    filtered_df['Mean_RPKM'] = filtered_df[['RPKM_EEF_54h_A', 'RPKM_EEF_54h_B']].mean(axis=1)

    # Create a new DataFrame with relevant columns
    result_df = filtered_df[['Unnamed: 0', 'Product Description', 'EEF_54h_A', 'EEF_54h_B', 'RPKM_EEF_54h_A', 'RPKM_EEF_54h_B','Mean_RPKM']]

    # Display the resulting DataFrame
    # Rename 'Unnamed: 0' to 'Gene ID' for clarity
    result_df.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)

    # print(result_df)
    return result_df

# getLiverStage()
