import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import chi2

colors = ['tab:green', 'tab:olive', 'tab:purple', 'tab:brown', 'tab:grey', 'black']
seq_names = ['Homo','N-term', 'C-term', 'M-Block', 'Cross']

parent_dir = os.path.join(os.getcwd(), 'RawData')


# Initialize dictionaries to hold traces and nucleation times
traces_0 = {i: [] for i in range(1, 6)} #nuleaction traces when no mRNA
traces_1 = {i: [] for i in range(1, 6)}
t_nuc_0 = {i: [] for i in range(1, 6)}
t_nuc_1 = {i: [] for i in range(1, 6)}


os.makedirs("t_nuc_values", exist_ok=True)
os.makedirs("time_traces", exist_ok=True)

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
            if NPro==3000:
                try:
                    df = pd.read_csv(os.path.join(folder_path, 'Largest_Cluster_Normalized.csv'),nrows=6000)
                    t = np.linspace(0, 1, len(df))  # normalized time
                    cluster_frac = df['ClusterFrac_Pro']
                    if NmRNA!= 0 :
                        mRNA_frac = df['ClusterFrac_mRNA']

                    # Determine t_nuc where ClusterFrac_Pro > 0.05
                    try:
                        t_nuc = t[np.where(cluster_frac > 0.05)[0][0]]
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
            trace_df.to_csv(f"time_traces/Seq{seq}_NmRNA{NmRNA}.csv", index=False)

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

from scipy.stats import expon
# Now plot
t_nuc_summary = []
for seq in range(1, 6):
    for NmRNA, traces, t_nuc_dict in [(0, traces_0, t_nuc_0), (1, traces_1, t_nuc_1)]:
        plt.figure(figsize=(6, 4))
        # Plot mean t_nuc with std
        t_nucs = np.array(t_nuc_dict[seq])

        # Set observation limit (censoring threshold), e.g. max simulation time
        T_obs = 1.0

        # Remove NaNs (uncensored)
        t_observed = t_nucs[~np.isnan(t_nucs)]
        if NmRNA ==0:
            l1= t_observed
        if NmRNA ==1:
            l2= t_observed
        if len(t_observed) > 0:
            k = len(t_observed)
            n = len(t_nucs)

            # Apply MLE formula
            lambda_hat = k / (np.sum(t_observed) + (n - k) * T_obs)
            mean = 1 / lambda_hat

            print(f"Estimated mean nucleation time (MLE): {mean:.4f}")

            # Sort t_observed times
            t_observed_sorted = np.sort(t_observed)
            ecdf = np.arange(1, k+1) / k

            plt.figure(figsize=(5,4))
            # Plot ECDF
            plt.step(t_observed_sorted, ecdf, where='post', label='Empirical CDF')

            # Plot fitted exponential CDF
            x = np.linspace(0, T_obs, 200)
            cdf_fit = expon.cdf(x, scale=mean)
            plt.plot(x, cdf_fit, label=f'Exponential Fit (mean={mean:.2f})')

            plt.xlabel('Simulation Time (normalized)', fontsize=13)
            plt.ylabel('Cumulative distribution function', fontsize=13)
            # plt.legend(title=f'{seq_names[seq-1]}')
            if NmRNA==1:
                plt.legend(title=f'{seq_names[seq-1]} with one polysome', fontsize=12)
            else:
                plt.legend(title=f'{seq_names[seq-1]} with no polysomes', fontsize=12)
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.tight_layout()
            plt.savefig(f't_nucs_Seq{seq}_NmRNA{NmRNA}.svg', format='svg',bbox_inches='tight', dpi=150)
            plt.close()

            SE = mean/np.sqrt(n)
            plt.figure(figsize=(5,4))
            # plt.axvline(mean,linewidth=0.5, color='black', linestyle='--', label=f't_nuc = {mean:.2f} ± {SE:.2f}')
            plt.axvspan(mean - SE, mean + SE, color='gray', alpha=0.3)

            # Save all t_nuc values per line
            df_t_nuc = pd.DataFrame({'t_nuc': t_observed_sorted})
            df_t_nuc.to_csv(f"t_nuc_values/Seq{seq}_NmRNA{NmRNA}.csv", index=False)

            # Save summary row
            t_nuc_summary.append({'Seq': seq, 'NmRNA': NmRNA, 'Mean_t_nuc': mean, 'SE_t_nuc': SE})

            if NmRNA==0:
            # Plot each replicate in transparent color
                for t, frac in traces[seq]:
                    plt.plot(t, frac, color='tab:blue',linewidth=0.5)
                    #plt.scatter(t, frac, color='tab:red', s=0.1, label='No Polysome' if 'No Polysome' not in plt.gca().get_legend_handles_labels()[1] else "")

                    custom_legend = [#Line2D([0], [0], color='black',linewidth=0, label=f'{seq_names[seq-1]}'),
                    #Line2D([0], [0], color='tab:blue', label='No polysome'),
                    Line2D([0], [0], linestyle='--', color='black',
                           label=fr'$t_{{\mathrm{{nuc}}}}$ = {mean:.2f} ± {SE:.2f}',
                           marker=None, linewidth=1.5)]
                    plt.legend(handles=custom_legend, loc='upper right', fontsize=12)

            # Example plotting logic
            if NmRNA == 1:
                
                # Combined data
                l_all = np.concatenate([l1, l2])

                # Compute MLE estimates
                λ1 = 1 / np.mean(l1)
                λ2 = 1 / np.mean(l2)
                λ_all = 1 / np.mean(l_all)

                # Compute log-likelihoods
                logL1 = len(l1) * np.log(λ1) - λ1 * np.sum(l1)
                logL2 = len(l2) * np.log(λ2) - λ2 * np.sum(l2)
                logL_all = len(l_all) * np.log(λ_all) - λ_all * np.sum(l_all)

                # LRT statistic
                Lambda = 2 * ((logL1 + logL2) - logL_all)

                # p-value
                p_value = 1 - chi2.cdf(Lambda, df=1)

                print(f"Seq: {seq_names[seq-1]}, Likelyhood ratio = {Lambda:.2f}, p = {p_value:.3f}")


                for t, frac, mRNA_frac in traces[seq]:
                    mask_high = mRNA_frac > 0.5
                    mask_low = ~mask_high

                    plt.plot(t, frac, color='tab:blue',linewidth=0.5)
                    # Plot high mRNA localization 
                    plot_masked_lines(t, frac, mask_high, color='tab:red', linewidth=0.5)

                    # Plot low mRNA localization 
                    #plot_masked_lines(t, frac, mask_low, color='tab:blue', linewidth=0.5)

                # Legend
                custom_legend = [#Line2D([0], [0], color='black',linewidth=0, label=f'{seq_names[seq-1]}'),
                   # Line2D([0], [0], color='tab:red', label='Polysome co-localized', linewidth=1.5),
                    #Line2D([0], [0], color='tab:blue', label='Polysome not co-localized', linewidth=1.5),
                    Line2D([0], [0], linestyle='--', color='black',
                           label=fr'$t_{{\mathrm{{nuc}}}}$ = {mean:.2f} ± {SE:.2f}, $p$ = {p_value:.2f}', linewidth=1.0)
                ]
                plt.legend(handles=custom_legend, loc='upper right', fontsize=13)

            # Formatting
            #plt.title(f'Seq {seq}, NmRNA={NmRNA}')
            plt.axvline(mean,linewidth=1.5, color='black', linestyle='--', label=f't_nuc = {mean:.2f} ± {SE:.2f}')
            plt.xlabel('Simulation Progress (normalized)', fontsize=15)
            plt.ylabel('Fraction of Largest Cluster', fontsize=15)
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            #plt.legend()
            plt.tight_layout()
            plt.savefig(f'Seq{seq}_NmRNA{NmRNA}.svg', format='svg',bbox_inches='tight', dpi=150)
            plt.close()
pd.DataFrame(t_nuc_summary).to_csv("t_nuc_summary.csv", index=False)
