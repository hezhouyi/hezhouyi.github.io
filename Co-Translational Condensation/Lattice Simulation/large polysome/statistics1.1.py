import os
import pandas as pd

# Parent directory containing your folders
parent_dir = os.getcwd()
parent_dir = os.path.join(parent_dir, 'RawData')

# Function to extract numbers from folder name
def extract_numbers(folder_name):
    parts = folder_name.split('_')  # Split the name by '_'
    
    try:
        seq = int(parts[1][3:])  # Extract the number from "seq1" (skip "seq")
        NmRNA = int(parts[2])    # Convert '5' to integer
        NPro = int(parts[3])     # Convert '1500' to integer
        replicate = int(parts[4])  # Convert '1' to integer
        #now I have many replicates, say 3, I want to read the same replicates as one, and then plot them later in a uncertainty way
        return seq, NmRNA, NPro, replicate
    except (IndexError, ValueError):
        print(f"Error parsing folder name: {folder_name}")
        return None

# Dictionary to store DataFrames for each unique seq number
PD_dict = {i: [] for i in range(1, 6)}  # Assuming seq numbers are 1-5
Drop_dict = {i: [] for i in range(1, 6)}  # Assuming seq numbers are 1-5

# Iterate through folders and extract numbers
for folder_name in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, folder_name)
    if os.path.isdir(folder_path) and folder_name.startswith("phase_seq"):
        result = extract_numbers(folder_name)
        if result:
            seq, NmRNA, NPro, replicate = result
            print(f"Folder: {folder_name}")
            print(f"  seq: {seq}, NmRNA: {NmRNA}, NPro: {NPro}, replicate: {replicate}")


            # Read only the third row of Phase_Diagram.csv into DataFrame
            phase_diagram = pd.read_csv(os.path.join(folder_path, 'Phase_Diagram.csv'), engine='python') 
            phase_diagram = phase_diagram.iloc[1:2]
            phase_diagram["replicate"] = replicate

            # Store the first row in the corresponding seq key in the dictionary
            if seq in PD_dict:
                PD_dict[seq].append(phase_diagram)

            droplet_info = pd.read_csv(os.path.join(folder_path, 'Droplet_info.csv'), engine='python',nrows=1) 

            # Add two more columns for NmRNA and NPro
            droplet_info["NmRNA"] = NmRNA
            droplet_info["NPro"] = NPro

            # Add a new column: N_mRNA / NmRNA
            if droplet_info["N_mRNA"].item() < 1:
                droplet_info["N_mRNA/NmRNA"] = 0
            else:
                droplet_info["N_mRNA/NmRNA"] = droplet_info["N_mRNA"] / NmRNA
            droplet_info["replicate"] = replicate

            # Store the first row in the corresponding seq key in the dictionary
            if seq in Drop_dict:
                Drop_dict[seq].append(droplet_info)

# After the loop, you will have 5 DataFrames (one for each seq number)

# To combine the data into single DataFrames for each seq number:
for seq, dfs in PD_dict.items():
    # Concatenate the list of DataFrames into one DataFrame
    PD_dict[seq] = pd.concat(dfs, ignore_index=True)

# Now, PD_dict contains 5 DataFrames for each seq (1-5). You can save them or process them further.
# For example, saving the DataFrames to CSV:
for seq, df in PD_dict.items():
    df=df.sort_values(by=[df.columns[1], df.columns[2]], ascending=[True, True])
    df.to_csv(f"phase_seq_{seq}.csv", index=False)


# To combine the data into single DataFrames for each seq number:
for seq, dfs in Drop_dict.items():
    # Concatenate the list of DataFrames into one DataFrame
    Drop_dict[seq] = pd.concat(dfs, ignore_index=True)

# Now, Drop_dict contains 5 DataFrames for each seq (1-5). You can save them or process them further.
# For example, saving the DataFrames to CSV:
for seq, df in Drop_dict.items():
    df=df.sort_values(by=[df.columns[-3], df.columns[-2]], ascending=[True, True])
    df.to_csv(f"droplet_seq_{seq}.csv", index=False)