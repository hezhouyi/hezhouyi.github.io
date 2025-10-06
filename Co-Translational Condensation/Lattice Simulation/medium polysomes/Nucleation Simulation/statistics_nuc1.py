import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

parent_dir = os.path.join(os.getcwd(), 'RawData')

# Initialize dictionaries to hold traces and nucleation times
traces_0 = {i: [] for i in range(1, 6)} #nuleaction traces when no mRNA
traces_1 = {i: [] for i in range(1, 6)}
t_nuc_0 = {i: [] for i in range(1, 6)}
t_nuc_1 = {i: [] for i in range(1, 6)}

# Function to extract folder parameters
def extract_numbers(folder_name):
    parts = folder_name.split('_')
    try:
        seq = int(parts[1][3:])
        NmRNA = int(parts[2])
        NPro = int(parts[3])
        replicate = int(parts[4])
        return seq, NmRNA, NPro, replicate
    except (IndexError, ValueError):
        return None

# Loop through folders and gather data
for folder_name in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, folder_name)
    if os.path.isdir(folder_path) and folder_name.startswith("phase_seq"):
        result = extract_numbers(folder_name)
        if result:
            seq, NmRNA, NPro, replicate = result
            try:
                df = pd.read_csv(os.path.join(folder_path, 'Largest_Cluster_Normalized.csv'), nrows=800)
                t = np.linspace(0, 1, len(df))  # normalized time
                cluster_frac = df['ClusterFrac_Pro']
                if NmRNA!= 0 :
                    mRNA_frac = df['ClusterFrac_mRNA']

                # Determine t_nuc where ClusterFrac_Pro > 0.2
                try:
                    t_nuc = t[np.where(cluster_frac > 0.2)[0][0]]
                except IndexError:
                    t_nuc = np.nan  # no crossing

                if NmRNA == 0:
                    traces_0[seq].append((t, cluster_frac))
                    t_nuc_0[seq].append(t_nuc)
                elif NmRNA == 1:
                    traces_1[seq].append((t, cluster_frac, mRNA_frac))
                    t_nuc_1[seq].append(t_nuc)

            except Exception as e:
                print(f"Error in {folder_name}: {e}")

# Save all traces (one file per seq and NmRNA) with shared time column
for NmRNA, traces_dict in [(0, traces_0), (1, traces_1)]:
    for seq, trace_list in traces_dict.items():
        if trace_list:
            # Assume all time arrays are the same length
            t = trace_list[0][0]  # shared time axis
            data = {'Time': t}

            # Add each trace as a new column
            if NmRNA==0:
                for idx, (time_arr, cluster_frac) in enumerate(trace_list):
                    data[f'Trace_{idx+1}'] = cluster_frac
                
            if NmRNA==1:
                for idx, (time_arr, cluster_frac, mRNA_frac) in enumerate(trace_list):
                    data[f'Trace_{idx+1}'] = cluster_frac
                    data[f'mRNATrace_{idx+1}'] = mRNA_frac

            trace_df = pd.DataFrame(data)
            trace_df.to_csv(f"Seq{seq}_NmRNA{NmRNA}_traces.csv", index=False)




# Save t_nuc summary and values
t_nuc_summary = []

for NmRNA, t_nuc_dict in [(0, t_nuc_0), (1, t_nuc_1)]:
    for seq, t_list in t_nuc_dict.items():
        t_arr = np.array(t_list)
        t_arr = t_arr[~np.isnan(t_arr)]
        if len(t_arr) > 0:
            mean = np.mean(t_arr)
            std = np.std(t_arr)
        else:
            mean = std = np.nan

        # Save all t_nuc values per line
        df_all = pd.DataFrame({'t_nuc': t_arr})
        df_all['Seq'] = seq
        df_all['NmRNA'] = NmRNA
        df_all.to_csv(f"Seq{seq}_NmRNA{NmRNA}_t_nuc_values.csv", index=False)

        # Save summary row
        t_nuc_summary.append({'Seq': seq, 'NmRNA': NmRNA, 'Mean_t_nuc': mean, 'Std_t_nuc': std})

# Save summary table
pd.DataFrame(t_nuc_summary).to_csv("t_nuc_summary.csv", index=False)

def plot_masked_lines(t, y, mask, color, **kwargs):
    """Plot a line but only where mask is True."""
    t = np.array(t)
    y = np.array(y)
    mask = np.array(mask)

    # Find continuous segments where mask is True
    i = 0
    while i < len(t):
        if mask[i]:
            j = i
            while j < len(t) and mask[j]:
                j += 1
            plt.plot(t[i:j], y[i:j], color=color, **kwargs)
            i = j
        else:
            i += 1

# Now plot
for seq in range(1, 6):
    for NmRNA, traces, t_nuc_dict in [(0, traces_0, t_nuc_0), (1, traces_1, t_nuc_1)]:
        plt.figure(figsize=(6, 4))
        # Plot mean t_nuc with std
        t_nucs = np.array(t_nuc_dict[seq])
        t_nucs = t_nucs[~np.isnan(t_nucs)]  # remove NaNs
        # if len(t_nucs) > 0:
        #     mean = np.mean(t_nucs)
        #     std = np.std(t_nucs)
        #     plt.axvline(mean, color='black', linestyle='--', label=f't_nuc = {mean:.2f} ± {std:.2f}')
        if len(t_nucs) > 0:
            mean = np.mean(t_nucs)
            std = np.std(t_nucs)
            plt.axvline(mean,linewidth=0.5, color='black', linestyle='--', label=f't_nuc = {mean:.2f} ± {std:.2f}')
            plt.axvspan(mean - std, mean + std, color='gray', alpha=0.3)


        if NmRNA==0:
        # Plot each replicate in transparent color
            for t, frac in traces[seq]:
                plt.plot(t, 100*frac, color='tab:blue',linewidth=0.5)
                #plt.scatter(t, 100 * frac, color='tab:red', s=0.1, label='No Polysome' if 'No Polysome' not in plt.gca().get_legend_handles_labels()[1] else "")

                custom_legend = [
                Line2D([0], [0], color='tab:blue', label='No polysome'),
                Line2D([0], [0], linestyle='--', color='black',
                       label=fr'$t_{{\mathrm{{nuc}}}}$ = {mean:.2f} ± {std:.2f}',
                       marker=None, linewidth=1.5)]
                plt.legend(handles=custom_legend, loc='upper right', fontsize=13)

        # if NmRNA==1:
        #     for t, frac, mRNA_frac in traces[seq]:
        #         #plt.plot(t, frac, color='tab:blue' if NmRNA == 0 else 'tab:red', alpha=0.3,linewidth=0.3)
        #         mask_high = mRNA_frac > 0.5 #mRNA is in the largest cluster
        #         mask_low = ~mask_high

        #         # Plot high values (blue): polysomes localized on condensates
        #         plt.scatter(t[mask_high], 100 * frac[mask_high], color='tab:blue', s=0.1)

        #         # Plot low values (red): polysomes not localized
        #         plt.scatter(t[mask_low], 100 * frac[mask_low], color='tab:red', s=0.1)

        #         custom_legend = [
        #         Line2D([0], [0], marker='o', color='w', label='Localized polysomes',
        #                markerfacecolor='tab:blue', markersize=5),
        #         Line2D([0], [0], marker='o', color='w', label='Non-localized polysomes',
        #                markerfacecolor='tab:red', markersize=5),
        #         Line2D([0], [0], linestyle='--', color='black',
        #                label=fr'$t_{{\mathrm{{nuc}}}}$ = {mean:.2f} ± {std:.2f}',
        #                marker=None, linewidth=1.5)
        #         ]

        #         plt.legend(handles=custom_legend, loc='upper right', fontsize=13)



        # Example plotting logic
        if NmRNA == 1:
            for t, frac, mRNA_frac in traces[seq]:
                mask_high = mRNA_frac > 0.5
                mask_low = ~mask_high

                plt.plot(t, 100*frac, color='tab:blue',linewidth=0.5)
                # Plot high mRNA localization 
                plot_masked_lines(t, 100 * frac, mask_high, color='tab:red', linewidth=0.5)

                # Plot low mRNA localization 
                #plot_masked_lines(t, 100 * frac, mask_low, color='tab:blue', linewidth=0.5)

            # Legend
            custom_legend = [
                Line2D([0], [0], color='tab:red', label='Polysome co-localized', linewidth=1.5),
                Line2D([0], [0], color='tab:blue', label='Polysome not co-localized', linewidth=1.5),
                Line2D([0], [0], linestyle='--', color='black',
                       label=fr'$t_{{\mathrm{{nuc}}}}$ = {mean:.2f} ± {std:.2f}', linewidth=1.0)
            ]
            plt.legend(handles=custom_legend, loc='upper right', fontsize=13)

        # Formatting
        #plt.title(f'Seq {seq}, NmRNA={NmRNA}')
        plt.xlabel('Simulation Progress (normalized)')
        plt.ylabel('Largest Cluster Percent (%)')
        plt.xlim(0, 1)
        plt.ylim(0, 100)
        #plt.legend()
        plt.tight_layout()
        plt.savefig(f'Seq{seq}_NmRNA{NmRNA}.svg', format='svg', dpi=150)
        plt.close()
