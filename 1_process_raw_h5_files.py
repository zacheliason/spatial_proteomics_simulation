import traceback

import matplotlib.pyplot as plt
import scipy.sparse as sp
import pandas as pd
import numpy as np
import subprocess
import shutil
import h5py
import json
import os
import re


root = '/Users/zacheliason/Documents/Work/payne/data_hallmarks'
pathologist_annotation_root = '/Users/zacheliason/Documents/Work/payne/20304456'

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

# Process files for each tissue
for tissue_id in ids:
    # Access spot annotations
    tissue_number = tissue_id

    # TODO: annotations
    # annotations_csv_path = os.path.join(pathologist_annotation_root, f'sample_{tissue_number}_pathologist_annotations.csv')
    # annotations_df = pd.read_csv(annotations_csv_path)

    # Process/save the spot positions list
    spot_positions_csv_path = os.path.join(root, tissue_id, "outs", "spatial", f"tissue_positions_list.csv")
    spot_positions = pd.read_csv(spot_positions_csv_path, header=None)

    # These headers come from documentation in https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    spot_positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']

    # Add annotations
    # spot_positions = spot_positions.merge(annotations_df, on='barcode', how='left')

    spot_positions.to_csv(spot_positions_csv_path, index=False)

    try:
        # Process/save the data inside each h5 file
        spatial_matrix_h5_path = os.path.join(root, tissue_id, "outs", f"filtered_feature_bc_matrix.h5")

        with h5py.File(spatial_matrix_h5_path, 'r') as h5_file:
            groups = list(h5_file.keys())

            # loop through each group and extract their datasets
            extracted_datasets = {}
            features = {}
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
                            if isinstance(dataset[feature], h5py._hl.dataset.Dataset):
                                extracted_datasets[feature] = dataset[feature][:]
                                features[feature] = dataset[feature][:]
                            else:
                                features[feature] = {}
                                for subfeature in dataset[feature]:
                                    extracted_datasets[subfeature] = dataset[feature][subfeature][:]
                                    features[feature][subfeature] = dataset[feature][subfeature][:]
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

    with open(os.path.join(root, tissue_id, f"h5_features.json"), 'w') as f:
        for feature in features:
            if type(features[feature]) == dict:
                for subfeature in features[feature]:
                    features[feature][subfeature] = [x.decode('utf-8') if isinstance(x, bytes) else x for x in features[feature][subfeature]]
            else:
                features[feature] = [x.decode('utf-8') if isinstance(x, bytes) else x for x in features[feature]]

        json.dump(features, f)

    # This creates a sparse matrix from the extracted datasets
    sparse_matrix = sp.csr_matrix((extracted_datasets['data'], extracted_datasets['indices'], extracted_datasets['indptr']), shape=extracted_datasets['shape'])

    matrix_indices = [x.decode('utf-8') for x in extracted_datasets['barcodes'].tolist()]
    matrix_columns = [x.decode('utf-8') for x in extracted_datasets['name'].tolist()]

    dense_df = pd.DataFrame(sparse_matrix.toarray())

    dense_df.columns = matrix_columns
    dense_df['barcode'] = matrix_indices

    dense_df.to_csv(os.path.join(root, tissue_id, "outs", f"processed_matrix.csv"), index=False)

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