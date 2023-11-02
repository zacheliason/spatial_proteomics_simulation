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

	df = df.rename(columns={'corners': 'grouping'})
	df['in_tissue'] = df.groupby('grouping')['in_tissue'].transform('max')
	df = df[df['grouping'] > 0]
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

	df = df.rename(columns={'corners': 'grouping'})
	df['in_tissue'] = df.groupby('grouping')['in_tissue'].transform('max')
	df = df[df['grouping'] > 0]
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

	df = df.rename(columns={'corners': 'grouping'})
	df['in_tissue'] = df.groupby('grouping')['in_tissue'].transform('max')
	df = df[df['grouping'] > 0]
	return df

root = '/Users/zacheliason/Documents/Work/payne/GSE208253_RAW'

# parameters
collapse_functions = [(collapse_spot_x_3, "collapse_x_3"), (collapse_spot_x_4, "collapse_x_4"), (collapse_spot_x_9, "collapse_x_9")]
perform_matrix_scaling = False
drop_incomplete_groupings = not perform_matrix_scaling

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]
for tissue_id in ids:
	spot_positions_csv_path = os.path.join(root, tissue_id, f"processed_{tissue_id}_spot_positions.csv")
	base_file_name = f"processed_{tissue_id}_matrix.csv"
	matrix_path = os.path.join(root, tissue_id, base_file_name)

	for collapse, suffix in collapse_functions:
		processed_matrix_path = os.path.join(root, tissue_id, base_file_name.replace(".csv", f"_{suffix}.csv"))

		loc_df = pd.read_csv(spot_positions_csv_path)
		loc_df = collapse(loc_df)

		# plt.scatter(loc_df['array_row'], loc_df['array_col'], s=5, c=loc_df['grouping'] % 5, cmap=plt.cm.coolwarm)
		# plt.show()

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
		matrix.to_csv(processed_matrix_path, index=False)
		print(f"Saved processed matrix to {processed_matrix_path}")
