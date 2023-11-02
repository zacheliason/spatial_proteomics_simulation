import matplotlib.pyplot as plt
import scipy.sparse as sp
import pandas as pd
import numpy as np
import subprocess
import shutil
import h5py
import os
import re


# Assign root as the RAW data directory (the first tar unzipped directory)
root = '/Users/zacheliason/Documents/Work/payne/GSE208253_RAW'

# If the directory has a space in it, replace it with an underscore
if " " in root:
    subprocess.run(['mv', root, root.replace(' ', '_')], check=True)
    root = root.replace(' ', '_')

# List files
filenames = [f for f in os.listdir(root) if os.path.isfile(os.path.join(root, f))]

# Unzip files if necessary
zipped_filenames = [f for f in filenames if f.endswith(".gz")]
if len(zipped_filenames) > 0:
    print("Unzipping files...")
    subprocess.run(["gunzip", f"{root}/*.gz"], check=True)

filenames = [f for f in filenames if " " not in f and not f.startswith(".")]
filenames = sorted(list(set(filenames)))

# Organizes files into directories based on the tissue_id if they haven't been already
if len(filenames) != 0:
    ids = []
    pattern = r"([\dA-Z]+_s\d+)_.*"
    for filename in filenames:
        match = re.match(pattern, filename)
        if match:
            id = match.group(1)
            ids.append(id)
            id_directory = os.path.join(root, id)

            if not os.path.exists(id_directory):
                os.mkdir(id_directory)

            shutil.move(os.path.join(root, filename), os.path.join(id_directory, filename))

# If files have already been organized, just get the ids from pre-existing directories
else:
    ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

# Process files for each tissue
for tissue_id in ids:
    # Process/save the spot positions list
    spot_positions_csv_path = os.path.join(root, tissue_id, f"{tissue_id}_tissue_positions_list.csv")
    spot_positions = pd.read_csv(spot_positions_csv_path, header=None)

    # These headers come from documentation in https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    spot_positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    spot_positions.to_csv(os.path.join(root, tissue_id, f"processed_{tissue_id}_spot_positions.csv"), index=False)

    try:
        # Process/save the data inside each h5 file
        spatial_matrix_h5_path = os.path.join(root, tissue_id, f"{tissue_id}_filtered_feature_bc_matrix.h5")

        with h5py.File(spatial_matrix_h5_path, 'r') as h5_file:
            groups = list(h5_file.keys())

            # loop through each group and extract their datasets
            extracted_datasets = {}
            for group_name in groups:
                group = h5_file[group_name]
                list_of_datasets = list(group.keys())

                # Loop through each dataset in the group
                for dataset_name in list_of_datasets:
                    dataset = group[dataset_name]

                    # Features are a special case because they are organized differently.
                    # This object stores 'genome', 'name', 'feature_type', and 'id'
                    # Of these, 'genome' and 'feature_type' are the same for all columns
                    # We will use the 'name' feature later as the columns of our df
                    if dataset_name == "features":
                        for feature in dataset:
                            extracted_datasets[feature] = dataset[feature][:]
                        continue

                    data = dataset[:]
                    extracted_datasets[dataset_name] = data

                    if dataset_name == "shape":
                        # Ensure that the shape is formatted correctly for scipy to read in our sparse matrix
                        if extracted_datasets['barcodes'].shape != data[0] and extracted_datasets['barcodes'].shape == data[1]:
                            extracted_datasets[dataset_name] = np.array([data[1], data[0]])
                        elif extracted_datasets['barcodes'].shape == data[0] and extracted_datasets['barcodes'].shape != data[1]:
                            extracted_datasets[dataset_name] = data
                        else:
                            # Haven't run into this case yet thank goodness
                            raise Exception(f"Error: {dataset_name} does not match barcodes shape")
    except Exception as e:
        print(f"Error: {e}, {dataset_name}")

    # This creates a sparse matrix from the extracted datasets
    sparse_matrix = sp.csr_matrix((extracted_datasets['data'], extracted_datasets['indices'], extracted_datasets['indptr']), shape=extracted_datasets['shape'])

    matrix_indices = [x.decode('utf-8') for x in extracted_datasets['barcodes'].tolist()]
    matrix_columns = [x.decode('utf-8') for x in extracted_datasets['name'].tolist()]

    dense_df = pd.DataFrame(sparse_matrix.toarray())

    dense_df.columns = matrix_columns
    dense_df['barcode'] = matrix_indices
    dense_df = dense_df[['barcode'] + matrix_columns[:-1]]

    dense_df.to_csv(os.path.join(root, tissue_id, f"processed_{tissue_id}_matrix.csv"), index=False)

    dense_df.index = dense_df['barcode']
    dense_df = dense_df.drop(columns=['barcode'])

    transcripts_per_spot = dense_df.sum(axis=1)
    transcripts_per_gene = dense_df.sum(axis=0)

    # count of non-zero genes across all positions
    unique_transcripts_per_spot = (dense_df > 0).sum(axis=1)

    # Print metrics for tissue
    print(f"Metrics for {tissue_id}:")
    print(f"  Total number of spots: {len(transcripts_per_spot)}")
    print(f"  Total number of genes: {len(transcripts_per_gene[transcripts_per_gene > 0])}")
    print(f"  Total number of transcripts: {transcripts_per_gene.sum()}")
    print(f"  Average number of transcripts per spot: {transcripts_per_spot.mean()}")
    print(f"  Average number of transcripts per gene: {transcripts_per_gene.mean()}")
    print(f"  Average number of unique transcripts per spot: {unique_transcripts_per_spot.mean()}")
    print()

    # Plot histograms
    # for i, r in dense_df.iterrows():
    #     r.hist(bins=100)
    #     plt.show()

    # unique_transcripts_per_spot.hist(bins=100)
    # plt.show()