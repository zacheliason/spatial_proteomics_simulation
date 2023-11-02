from matplotlib import pyplot as plt
from itertools import count
import pandas as pd
import numpy as np
import os


# Collapse spots into larger groupings of 9 spots
def collapse_spot_x_9(df): # large squares
	def identify_corners(row, counter):
		x = row['array_row']
		y = row['array_col']

		if y % 3 == 0 and x % 3 == 0:
			return next(counter)
		else:
			return 0

	id_counter = count(1)
	df['corners'] = df.apply(identify_corners, args=(id_counter,), axis=1)

	corner_indices = df[df['corners'] > 0].index.values.tolist()
	for i, r in df.iterrows():
		if i not in corner_indices:
			continue

		x = r['array_row']
		y = r['array_col']
		c = r['corners']

		df.loc[(df['array_row'] == x-2) & (df['array_col'] == y), 'corners'] = c
		df.loc[(df['array_row'] == x-1) & (df['array_col'] == y-1), 'corners'] = c
		df.loc[(df['array_row'] == x-1) & (df['array_col'] == y+1), 'corners'] = c
		df.loc[(df['array_row'] == x) & (df['array_col'] == y+2), 'corners'] = c
		df.loc[(df['array_row'] == x) & (df['array_col'] == y-2), 'corners'] = c
		df.loc[(df['array_row'] == x+1) & (df['array_col'] == y+1), 'corners'] = c
		df.loc[(df['array_row'] == x+1) & (df['array_col'] == y-1), 'corners'] = c
		df.loc[(df['array_row'] == x+2) & (df['array_col'] == y), 'corners'] = c

	# df = df.rename(columns={'corners': 'grouping'})
	df['grouping'] = df['corners'].astype(str) + "x9"

	df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('max')
	df = df[df['corners'] > 0]
	return df


# Collapse spots into larger groupings of 4 spots
def collapse_spot_x_4(df): # small squares
	def identify_corners(row, counter):
		x = row['array_row']
		y = row['array_col']

		if y % 4 == 2 and x % 4 == 0:
				return next(counter)
		elif y % 4 == 0 and x % 4 == 2:
				return next(counter)
		else:
			return 0

	id_counter = count(1)
	df['corners'] = df.apply(identify_corners, args=(id_counter,), axis=1)

	corner_indices = df[df['corners'] > 0].index.values.tolist()
	for i, r in df.iterrows():
		if i not in corner_indices:
			continue

		x = r['array_row']
		y = r['array_col']
		c = r['corners']

		df.loc[(df['array_row'] == x) & (df['array_col'] == y - 2), 'corners'] = c
		df.loc[(df['array_row'] == x + 1) & (df['array_col'] == y - 1), 'corners'] = c
		df.loc[(df['array_row'] == x - 1) & (df['array_col'] == y - 1), 'corners'] = c

	# df = df.rename(columns={'corners': 'grouping'})
	df['grouping'] = df['corners'].astype(str) + "x4"

	df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('max')
	df = df[df['corners'] > 0]
	return df


# Collapse spots into larger groupings of 3 spots
def collapse_spot_x_3(df): # triangles
	def identify_corners(row, counter):
		x = row['array_row']
		y = row['array_col']

		if y % 4 == 1:
			if x % 2 == 1 and x % 3 != 2:
				return next(counter)
			else:
				return 0
		elif y % 4 == 3:
			if x % 2 == 1 and x % 3 != 1:
				return next(counter)
			else:
				return 0
		else:
			return 0

	id_counter = count(1)
	df['corners'] = df.apply(identify_corners, args=(id_counter,), axis=1)

	corner_indices = df[df['corners'] > 0].index.values.tolist()
	for i, r in df.iterrows():
		if i not in corner_indices:
			continue

		x = r['array_row']
		y = r['array_col']
		c = r['corners']

		if y % 4 == 3:
			if x % 3 == 0:  # corner is right of triangle
				df.loc[(df['array_row'] == x - 2) & (df['array_col'] == y), 'corners'] = c
				df.loc[(df['array_row'] == x - 1) & (df['array_col'] == y - 1), 'corners'] = c
			else:           # corner is top of triangle
				df.loc[(df['array_row'] == x + 1) & (df['array_col'] == y - 1), 'corners'] = c
				df.loc[(df['array_row'] == x - 1) & (df['array_col'] == y - 1), 'corners'] = c
		else:
			if x % 3 == 0:  # corner is left of triangle
				df.loc[(df['array_row'] == x + 2) & (df['array_col'] == y), 'corners'] = c
				df.loc[(df['array_row'] == x + 1) & (df['array_col'] == y - 1), 'corners'] = c
			else:           # corner is top of triangle
				df.loc[(df['array_row'] == x + 1) & (df['array_col'] == y - 1), 'corners'] = c
				df.loc[(df['array_row'] == x - 1) & (df['array_col'] == y - 1), 'corners'] = c

	# df = df.rename(columns={'corners': 'grouping'})
	df['grouping'] = df['corners'].astype(str) + "x3"

	df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('max')
	df = df[df['corners'] > 0]
	return df

root = '/Users/zacheliason/Documents/Work/payne/GSE208253_RAW'

# parameters
collapse_functions = [(collapse_spot_x_3, 3), (collapse_spot_x_4, 4), (collapse_spot_x_9, 9)]
perform_matrix_scaling = False
drop_incomplete_groupings = not perform_matrix_scaling

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]
for tissue_id in ids:
	spot_positions_csv_path = os.path.join(root, tissue_id, f"processed_{tissue_id}_spot_positions.csv")
	base_file_name = f"processed_{tissue_id}_matrix.csv"
	matrix_path = os.path.join(root, tissue_id, base_file_name)

	spot_positions_df = pd.read_csv(spot_positions_csv_path)

	for collapse, grouping_factor in collapse_functions:
		simulated_matrix_path = os.path.join(root, tissue_id, 'simulated_matrices', base_file_name.replace(".csv", f"_collapsed_x_{grouping_factor}.csv").replace("processed_", ""))
		simulated_loc_path = os.path.join(root, tissue_id, 'simulated_matrices', f"{tissue_id}_spot_positions_collapsed_x_{grouping_factor}.csv")
		if not os.path.exists(os.path.dirname(simulated_matrix_path)):
			os.makedirs(os.path.dirname(simulated_matrix_path))

		loc_df = spot_positions_df.copy()
		loc_df = collapse(loc_df)
		# create new_x which is average of array_row for each grouping
		loc_df['centroid_row'] = loc_df.groupby('grouping')['array_row'].transform('mean')
		loc_df['centroid_col'] = loc_df.groupby('grouping')['array_col'].transform('mean')
		loc_df = loc_df.sort_values(by=['centroid_row', 'centroid_col'])

		loc_df.to_csv(simulated_loc_path, index=False)
		loc_df = loc_df[loc_df['in_tissue'] == 1]

		# plot the groupings
		x_dim = 5
		y_dim = loc_df['array_col'].max() / loc_df['array_row'].max() * x_dim

		plt.figure(figsize=(x_dim, y_dim))
		plt.scatter(loc_df['array_row'], loc_df['array_col'], s=2, c=loc_df['corners'] % 8, cmap=plt.cm.coolwarm)
		plt.scatter(loc_df['centroid_row'], loc_df['centroid_col'], s=4, c="black", label="New Grouping Center")
		for i, r in loc_df.iterrows():
			plt.plot([r['centroid_row'], r['array_row']], [r['centroid_col'], r['array_col']], c='black', alpha=.5, linewidth=.7)

		plt.savefig(os.path.join(root, tissue_id, 'simulated_matrices', f"{tissue_id}_collapsed_x_{grouping_factor}.png"))
		plt.close()

		loc_df = loc_df[['barcode', 'grouping']]

		# Read in expression matrix
		df = pd.read_csv(matrix_path)

		# Merge expression matrix with new grouping ids
		df = df.merge(loc_df, on='barcode', how='left')

		# Count how many spots are in each grouping
		df['num_per_grouping'] = df.groupby('grouping')['grouping'].transform('count')

		# Drop any groupings which are missing spots in the tissue
		if drop_incomplete_groupings:
			df = df[df['num_per_grouping'] == max(df['num_per_grouping'])]

		# Save scaling factors for each grouping (based on how many spots each is missing)
		num_per_grouping_dict = df.drop_duplicates(subset='grouping').set_index('grouping')['num_per_grouping'].to_dict()
		scalar_dict = {k: max(num_per_grouping_dict.values()) / v for k, v in num_per_grouping_dict.items()}

		# Drop unnecessary columns and sum gene expression data by new grouping id
		df = df.drop(columns=['barcode', 'num_per_grouping'])
		matrix = df.groupby('grouping').sum()

		if perform_matrix_scaling:
			scaled_matrix = matrix.copy().astype(float)

			# create numpy float array from dictionary values sorted by key
			scalar_array = np.array(list(dict(sorted(scalar_dict.items())).values()), dtype=float)
			scalar_array = scalar_array[:, np.newaxis]

			# create numpy float array from dataframe values
			scaled_matrix_reset = scaled_matrix.reset_index(drop=True)
			scaled_array = scaled_matrix_reset.to_numpy()

			# scale each row by the scalar value for that grouping
			scaled_array *= scalar_array

			# recreate dataframe with scaled values
			scaled_matrix_reset = pd.DataFrame(data=scaled_array, columns=scaled_matrix.columns)
			scaled_matrix_reset.index = scaled_matrix.index
			matrix = scaled_matrix_reset

		# format index groupings for output
		matrix['groupings'] = matrix.index
		matrix = matrix.reset_index(drop=True)[['groupings'] + matrix.columns[:-1].tolist()]
		matrix.to_csv(simulated_matrix_path, index=False)
		print(f"Saved simulated x {grouping_factor} matrix to {simulated_matrix_path}")
