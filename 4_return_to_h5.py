import scipy.sparse as sp
import numpy as np
import pandas as pd
import h5py
import json
import os
import re


def matrix_to_h5(matrix_path, features_path, output_file_path):
    # Load the dense matrix and features json
    dense_df = pd.read_csv(matrix_path)
    with open(features_path, 'r') as f:
        features = json.load(f)

        # Convert all features to byte strings
        for feature in features:
            if type(features[feature]) == dict:
                for subfeature in features[feature]:
                    features[feature][subfeature] = [x.encode('utf-8') for x in features[feature][subfeature]]
            else:
                features[feature] = [x.encode('utf-8') for x in features[feature]]

    # Get the matrix indices and columns
    matrix_indices = dense_df['barcode'].values
    matrix_columns = [x for x in dense_df.columns.values if x != 'barcode']

    # Get the raw values of the matrix
    raw_values = dense_df.drop(columns=['barcode']).values

    # Convert matrix indices back to byte strings
    matrix_indices_bytes = [x.encode('utf-8') for x in matrix_indices]

    # Create a sparse matrix from the raw values
    matrix_data = sp.csr_matrix(raw_values)

    # Create an HDF5 file
    with h5py.File(output_file_path, 'w') as h5_file:
        # Create a group for the matrix data
        group = h5_file.create_group('matrix')

        # Store matrix data
        group.create_dataset('data', data=matrix_data.data)
        group.create_dataset('indices', data=matrix_data.indices)
        group.create_dataset('indptr', data=matrix_data.indptr)
        group.create_dataset('shape', data=matrix_data.shape)

        # Store barcodes and features
        group.create_dataset('barcodes', data=np.array(matrix_indices_bytes, dtype='S'))
        feature_group = group.create_group('features')
        for feature in features:
            if type(features[feature]) == dict:
                subgroup = feature_group.create_group(feature)
                for subfeature in features[feature]:
                    subgroup.create_dataset(subfeature, data=np.array(features[feature][subfeature], dtype='S'))
            else:
                feature_group.create_dataset(f'{feature}', data=np.array(features[feature], dtype='S'))

    print(f'Finished writing HDF5 file to {output_file_path}')
    os.remove(matrix_path)



root = '/Users/zacheliason/Documents/Work/payne/data_hallmarks'

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

# filter all ids by regex pattern
pattern = re.compile(r'^\w+\d*_x_\d')
ids = [f for f in ids if pattern.match(f)]

for tissue_id in ids:
    if tissue_id == "Lung7_x_4":
        continue


    tissue_dir = os.path.join(root, tissue_id)

    matrix_path = os.path.join(tissue_dir, 'outs', 'processed_matrix.csv')
    features_path = os.path.join(tissue_dir, "h5_features.json")

    reverted_h5_path = os.path.join(tissue_dir, 'outs', 'filtered_feature_bc_matrix.h5')

    if os.path.exists(reverted_h5_path):
        os.remove(reverted_h5_path)

    matrix_to_h5(matrix_path=matrix_path, features_path=features_path, output_file_path=reverted_h5_path)