from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
from scipy.spatial import ConvexHull
from itertools import count
import pandas as pd
import numpy as np
import shutil
import glob
import os


# Collapse spots into larger groupings of 9 spots
def collapse_spot_x_9(df, drop_incomplete_pairings=True): # large squares
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

	if drop_incomplete_pairings:
		df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('min')

	df = df[df['corners'] > 0]
	return df


# Collapse spots into larger groupings of 4 spots
def collapse_spot_x_4(df, drop_incomplete_pairings=True): # small squares
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

	if drop_incomplete_pairings:
		df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('min')
	df = df[df['corners'] > 0]
	return df


# Collapse spots into larger groupings of 3 spots
def collapse_spot_x_3(df, drop_incomplete_pairings=True): # triangles
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

	if drop_incomplete_pairings:
		df['in_tissue'] = df.groupby('corners')['in_tissue'].transform('min')
	df = df[df['corners'] > 0]
	return df

root = '/Users/zacheliason/Documents/Work/payne/data_hallmarks'

# collapse_functions = [(collapse_spot_x_3, 3), (collapse_spot_x_4, 4), (collapse_spot_x_9, 9)]
collapse_functions = [(collapse_spot_x_4, 4), (collapse_spot_x_9, 9)]
drop_incomplete_groupings = True

ids = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

for tissue_id in ids:
	for collapse, grouping_factor in collapse_functions:

		# copy entire tissue_id directory
		new_tissue_directory = os.path.join(root, tissue_id + "_x_" + str(grouping_factor))
		if os.path.exists(new_tissue_directory):
			shutil.rmtree(new_tissue_directory)

		shutil.copytree(os.path.join(root, tissue_id), new_tissue_directory)

		spot_positions_csv_path = os.path.join(new_tissue_directory, "outs", "spatial", f"tissue_positions_list.csv")
		barcode_to_groupings_path = os.path.join(new_tissue_directory, f"barcodes_to_collapsed_spot_groupings.csv")
		matrix_path = os.path.join(new_tissue_directory, "outs", "processed_matrix.csv")

		spot_positions_df = pd.read_csv(spot_positions_csv_path)

		loc_df = spot_positions_df.copy()
		# drop any rows with strings in 'array_row' column (in case there are multiple column headers)

		# check if any string types in array_row column
		if loc_df['array_row'].dtype == np.object:
			loc_df = loc_df[~loc_df['array_row'].str.contains("[a-zA-Z]").fillna(False)]
		loc_df[['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']] = loc_df[['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']].astype(int)

		loc_df = collapse(loc_df)
		# create new_x which is average of array_row for each grouping
		loc_df['centroid_row'] = loc_df.groupby('grouping')['array_row'].transform('mean')
		loc_df['centroid_col'] = loc_df.groupby('grouping')['array_col'].transform('mean')
		loc_df['pxl_centroid_row'] = loc_df.groupby('grouping')['pxl_row_in_fullres'].transform('mean')
		loc_df['pxl_centroid_col'] = loc_df.groupby('grouping')['pxl_col_in_fullres'].transform('mean')
		loc_df = loc_df.sort_values(by=['centroid_row', 'centroid_col'])

		loc_df.to_csv(barcode_to_groupings_path, index=False)

		# This df is only used to save the new spot positions to a csv
		# This returns it to its original format
		save_loc_df = loc_df.drop(columns=['array_row', 'array_col', 'corners', 'barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'])
		save_loc_df = save_loc_df.rename(columns={'centroid_row': 'array_row', 'centroid_col': 'array_col', 'grouping': 'barcode', 'pxl_centroid_row': 'pxl_row_in_fullres', 'pxl_centroid_col': 'pxl_col_in_fullres'})
		save_loc_df = save_loc_df[['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]
		save_loc_df = save_loc_df.drop_duplicates()

		save_loc_df.to_csv(spot_positions_csv_path, index=False)

		loc_df = loc_df[loc_df['in_tissue'] == 1]

		# TODO: add annotation column to processed_spot_positions.csv
		# loc_df['pathologist_annotation'] = loc_df['pathologist_annotation'].fillna("No annotation")
		# loc_df['cluster_annotation'] = loc_df['pathologist_annotation'].fillna("None")

		# If figure doesn't already exist, plot the new spot groupings overlaid on the tissue
		figpath = os.path.join(new_tissue_directory, f"collapsed_spots_map.png")
		if not os.path.exists(figpath):
			# plot the groupings
			x_dim = 5
			y_dim = (loc_df['array_col'].max() / loc_df['array_row'].max() * x_dim) + 10

			plt.style.use('seaborn-whitegrid')
			plt.figure(figsize=(x_dim, y_dim))
			fig, ax = plt.subplots()

			plt.legend()
			for grouping in loc_df['grouping'].unique():
				grouping_df = loc_df[loc_df['grouping'] == grouping].reset_index(drop=True)
				if len(grouping_df) < 3:
					continue

				hull = ConvexHull(grouping_df[['array_row', 'array_col']])
				hull_vertices = grouping_df.iloc[hull.vertices][['array_row', 'array_col']].to_numpy()

				polygon = Polygon(xy=hull_vertices, closed=True, facecolor='black', edgecolor='black', alpha=.4)
				ax.add_patch(polygon)

			# TODO remove this line
			loc_df['pathologist_annotation'] = "No annotation"

			cmap = plt.get_cmap('tab10')
			for annotation, color in zip(loc_df['pathologist_annotation'].unique(), cmap.colors):
				if annotation == "No annotation":
					color = "grey"
				category_data = loc_df[loc_df['pathologist_annotation'] == annotation]
				plt.scatter(category_data['array_row'], category_data['array_col'], s=4, c=color, label=annotation)

			ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=False)

			plt.savefig(figpath, bbox_inches='tight', dpi=300)
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

		# # Save scaling factors for each grouping (based on how many spots each is missing)
		# num_per_grouping_dict = df.drop_duplicates(subset='grouping').set_index('grouping')['num_per_grouping'].to_dict()
		# scalar_dict = {k: max(num_per_grouping_dict.values()) / v for k, v in num_per_grouping_dict.items()}

		# Drop unnecessary columns and sum gene expression data by new grouping id
		df = df.drop(columns=['barcode', 'num_per_grouping'])

		# mask_df marks whether there is any expression for a given gene in a given spot (1 or 0)
		gene_columns = df.columns[df.columns != 'grouping']
		mask_df = df.copy()
		mask_df[gene_columns] = (mask_df[gene_columns] != 0).astype(int)

		# matrix is the expression matrix collapsed by new grouping ids
		matrix = df.groupby('grouping').sum()

		# found_in_n_spots marks how many spots a given gene is found in among each collapsed spot grouping
		found_in_n_spots = mask_df.groupby('grouping').sum()

		# Drop any genes which are found in DROPPING_THRESHOLD or fewer spots
		DROPPING_THRESHOLD = 2

		nonzero = matrix > 0
		unique_nonzero_cols_per_row = nonzero.sum(axis=1)

		fig = plt.figure(figsize=(10, 5))
		plt.style.use("ggplot")
		plt.hist(unique_nonzero_cols_per_row, bins=100)

		# plot horizontal line for mean and median
		plt.axvline(x=np.mean(unique_nonzero_cols_per_row), color='blue', linestyle='dashed', linewidth=1, label='mean')
		plt.axvline(x=np.median(unique_nonzero_cols_per_row), color='green', linestyle='dashed', linewidth=1, label='median')
		plt.legend()

		plt.title(f"{tissue_id} number of unique transcripts per collapsed spot\n before dropping counts found in {DROPPING_THRESHOLD} or fewer spots:")
		plt.ylabel("unique transcripts per collapsed spot")

		plt.savefig(os.path.join(new_tissue_directory, f"unique_transcripts_NO_DROP_x_{grouping_factor}.png"), bbox_inches='tight', dpi=300)
		plt.clf()
		plt.close()

		# Drop counts which are only found in DROPPING_THRESHOLD or fewer spots
		mask_matrix = found_in_n_spots.apply(lambda x: x < (grouping_factor // 4) + DROPPING_THRESHOLD)
		matrix = matrix.mask(mask_matrix)
		matrix = matrix.fillna(0)

		nonzero = matrix > 0
		unique_nonzero_cols_per_row = nonzero.sum(axis=1)

		fig = plt.figure(figsize=(10, 5))
		plt.style.use("ggplot")
		plt.hist(unique_nonzero_cols_per_row, bins=100)

		plt.axvline(x=np.mean(unique_nonzero_cols_per_row), color='blue', linestyle='dashed', linewidth=1, label='mean')
		plt.axvline(x=np.median(unique_nonzero_cols_per_row), color='green', linestyle='dashed', linewidth=1, label='median')
		plt.legend()

		plt.title(f"{tissue_id} number of unique transcripts per collapsed spot\n after dropping counts found in {DROPPING_THRESHOLD} or fewer spots:")
		plt.ylabel("unique transcripts per collapsed spot")

		plt.savefig(os.path.join(new_tissue_directory, f"unique_transcripts_WITH_DROP_x_{grouping_factor}.png"), bbox_inches='tight', dpi=300)

		sum_per_row = matrix.sum(axis=1)
		plt.hist(sum_per_row, bins=100)

		# format index groupings for output
		matrix['barcode'] = matrix.index
		matrix = matrix.reset_index(drop=True)
		matrix.to_csv(matrix_path, index=False)
		print(f"Saved simulated x {grouping_factor} matrix to {matrix_path}")