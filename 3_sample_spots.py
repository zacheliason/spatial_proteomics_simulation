import pandas as pd
import numpy as np
import os
import re


# spataial proteomics in the proteomics world is more rare than it is in the lipidomics world because of these
# pros and cons (talk about how with MALDI lipidomics they raster a single spot into a single spectrum
# whereas with MALDI proteomics they run a 1 hour test on each LCM section


#Q what to do about spots that fall out of the tissue but which are inside the grouping?
	# drop the grouping? scale the counts up to the grouping size?
	# need to figure out normalization for the matrices


def random_sample(df, n=20, filter_column=None):
	df = df[df['in_tissue'] == 1]
	if filter_column is None:
		random_indices = np.random.choice(df.index, size=n, replace=False)
	else:
		random_groupings = np.random.choice(df.drop_duplicates(subset=filter_column)[filter_column], size=n,
		                                    replace=False)
		random_indices = df[df[filter_column].isin(random_groupings)].index

	return df.loc[random_indices].sort_values(by=['array_row', 'array_col'])


root = '/Users/zacheliason/Documents/Work/payne/GSE208253_RAW'

# parameters
perform_matrix_scaling = False
drop_incomplete_groupings = not perform_matrix_scaling

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]
for tissue_id in ids:
	# Add collapsed spot position files
	spot_positions_files = sorted([os.path.join(root, tissue_id, 'simulated_matrices', f) for f in os.listdir(os.path.join(root, tissue_id, 'simulated_matrices')) if 'spot_positions' in f])
	# Add collapsed matrix files
	matrix_files = sorted([os.path.join(root, tissue_id, 'simulated_matrices', f) for f in os.listdir(os.path.join(root, tissue_id, 'simulated_matrices')) if 'matrix' in f])
	# Add original processed spot position file
	spot_positions_files.append(os.path.join(root, tissue_id, f"processed_{tissue_id}_spot_positions.csv"))
	# Add original processed matrix file
	matrix_files.append(os.path.join(root, tissue_id, f"processed_{tissue_id}_matrix.csv"))

	# Zip together spot positions and matrix files
	file_pairs = list(zip(spot_positions_files, matrix_files))

	# Iterate through each spot positions/matrix file pair and create new sampled matrix files from each
	for spot_positions_path, matrix_path in file_pairs:
		matrix_base_name = os.path.basename(matrix_path)
		filter_settings = []

		spot_positions_df = pd.read_csv(spot_positions_path)
		matrix = pd.read_csv(matrix_path)

		grouping_col = matrix.columns[0]

		# SAMPLE key positions list
		key_positions_list = [(x,y) for x in range(1, 9) for y in range(1, 9)]

		key_filtered_spot_positions_df = spot_positions_df[spot_positions_df.apply(lambda x: (x['array_row'], x['array_col']) in key_positions_list, axis=1)]
		key_filtered_matrix_pathname = os.path.join(root, tissue_id, 'simulated_matrices', matrix_base_name.replace(".csv", "_key_filtered.csv"))
		filter_settings.append(key_filtered_spot_positions_df, key_filtered_matrix_pathname)

		# Make random filtered sample
		random_filtered_spot_positions_df = random_sample(spot_positions_df, n=20, filter_column='grouping')
		random_filtered_matrix_pathname = os.path.join(root, tissue_id, 'simulated_matrices', matrix_base_name.replace(".csv", "_random_filtered.csv"))
		filter_settings.append(random_filtered_spot_positions_df, random_filtered_matrix_pathname)

		for filtered_positions_df, matrix_filepath in filter_settings:
			filtered_matrix = matrix[matrix[grouping_col].isin(filtered_positions_df[grouping_col])]
			filtered_matrix.to_csv(matrix_filepath, index=False)
