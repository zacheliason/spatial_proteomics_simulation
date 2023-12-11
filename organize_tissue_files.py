import subprocess
import shutil
import re
import os


# Assign root as the RAW data directory (the first tar unzipped directory)
root = '/Users/zacheliason/Documents/Work/payne/hallmarks/data_hallmarks'

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
