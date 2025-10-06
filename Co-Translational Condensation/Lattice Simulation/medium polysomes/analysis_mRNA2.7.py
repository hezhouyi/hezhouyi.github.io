#!/usr/bin/env python
# coding: utf-8
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def get_variable(variable_name,filename='y_runkey.txt'):
    values = []
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Flag to indicate whether to collect lines
    collect_lines = False
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace
        if collect_lines:
            values.extend([int(x) if x.isdigit() else x for x in line.split()])  # Split the line and extend the values list
            break  # Stop collecting lines when an empty line is encountered
        elif line == variable_name:
            collect_lines = True

    # If the list has only one element, return it as a single value
    if len(values) == 1:
        return values[0]
    # Otherwise, return the list of values
    else:
        return values
    
# Define the sigmoid function
def sigmoid(x, a, b, c):
    return a / (1 + np.exp(b * (x - c)))

def DouSig(x, a1, a2, b1, b2, c1, c2):
    return -a1 /(1 + np.exp(b1 * (x - c1)))+ a2 /(1 + np.exp(b2 * (x - c2)))


#Read Keywords:
BoxSize = get_variable('BoxSize')
if BoxSize[0]==BoxSize[2]:
    if BoxSize[0]==100:
        RD=87
    elif BoxSize[0]==200:
        RD=174

    elif  BoxSize[0]==400:
        RD=174
        print('Your boxsize 400, which is not ideal, but lets continue')
    else:
        print('Your boxsize in wrong!!!')
        sys.exit()

NkN = get_variable('NkN')
Nk = get_variable('Nk')
Vk = get_variable('Vk')
Vak = get_variable('Vak')

if NkN>1:
    Vk_m=max(Vk)
else:
    Vk_m=Vk
Tmax = get_variable('Tmax')  # strip any trailing newline characters
E_l = get_variable('E_LargestCluster')
m = np.dot(Vk, Nk)
n = int(Tmax/E_l + 1)
V_B=BoxSize[0]*BoxSize[1]*BoxSize[2]
StartDroplet=get_variable('StartDroplet')
StartSlab=get_variable('StartSlab')
C_µM = 10**7/(6.02*27)  # Extra 10^6 is to convert M to µM, one lattice site to 3nm
#plot Number of Clusters of all data in subfolders
x=[]
cluster_hist=[]
#plot Smallest Cluster Size
data=pd.read_csv('A_ClusterHist.txt',engine='python', sep='\s+', header=None, names=['Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6', 'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12', 'Column13', 'Column14', 'Column15', 'Column16', 'Column17', 'Column18', 'Column19', 'Column20', 'Column21'])
others = data.Column5+data.Column6+data.Column7+data.Column8+data.Column9+data.Column10+data.Column11+data.Column12+data.Column13+data.Column14+data.Column15+data.Column16+data.Column17+data.Column18+data.Column19+data.Column20+data.Column21
n_avg=(len(data.Column1)-1)//10 # use last 10% points
for i in range(21):
    cluster_hist.append(sum(data.iloc[-n_avg:, i]) / n_avg)
    x.append(str(i+1))
x[20]='21+'
plt.bar(x, cluster_hist)
for a,b in zip(x,cluster_hist):
    plt.text(a, b+0.05, '%.1f' % b, ha='center', va= 'bottom',fontsize=7)

plt.xlabel('Cluster Size')
plt.ylabel('Distribution')
plt.savefig('Cluster Size Distribution.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

plt.figure(figsize=(6.3, 3.6))
plt.plot(data.Column1,'g', label='Monomer ~ '+ '%.1f' % cluster_hist[0])
plt.plot(data.Column2, 'b', label='Dimer ~ '+ '%.1f' % cluster_hist[1])
plt.plot(data.Column3, 'y', label='Trimer ~ '+ '%.1f' % cluster_hist[2])
plt.plot(data.Column4, 'r', label='Tetramer ~ '+ '%.1f' % cluster_hist[3])
plt.plot(others ,'black', label='Others ~ '+ '%.1f' % sum(cluster_hist[4:])) 
plt.xlabel('Simulation Process')
plt.ylabel('Number of Clusters')
plt.legend(loc='upper right')
plt.savefig('Cluster_size_time_course.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

#plot Largest Cluster Size
data=pd.read_csv('A_LargestCluster.txt',engine='python',sep='\s+',  header=None)
n_avg=(len(data.iloc[:,0])-1)//2 #Use last 50% as average value
avg1 = sum(data.iloc[-n_avg:, 0]) / n_avg #python starts with index 0!
avg2 = sum(data.iloc[-n_avg:, NkN+1]) / n_avg
#avgCore = sum(data.iloc[-n_avg:, 3*(NkN+1)]) / n_avg 
#avgCore1 = sum(data.iloc[-n_avg:, 4*(NkN+1)]) / n_avg
DropProteinN= []
avg_pro1=sum(data.iloc[-n_avg:, 2]) / n_avg
avg_pls1=sum(data.iloc[-n_avg:, 1]) / n_avg
avg_pro2=sum(data.iloc[-n_avg:, NkN+1+2]) / n_avg
avg_pls2=sum(data.iloc[-n_avg:, NkN+1+1]) / n_avg
avg_pro3=sum(data.iloc[-n_avg:, 2*(NkN+1)+2]) / n_avg
avg_pls3=sum(data.iloc[-n_avg:, 2*(NkN+1)+1]) / n_avg
N_Pro=avg_pro1
if avg_pro2>5 or avg_pls2>2:
    N_Pro+=avg_pro2
    if avg_pro3>5 or avg_pls3>2:
        N_Pro+=avg_pro3
DropProteinN.append(avg_pls1)
DropProteinN.append(N_Pro)


plt.scatter(range(len(data.iloc[:, 1])),data.iloc[:,0],color='r', label='Largest Cluster ~ '+ '%.2f' % avg1,s=0.5)

#plt.plot(data.iloc[:,0],'r', linestyle=':',label='Largest Cluster ~ '+ '%.2f' % avg1)

#linear regression fitting to test if converged
from scipy.stats import linregress
x_values=range(n-n_avg+1,n+1)
y_values=data.iloc[-n_avg:, 0]
slope, intercept, r_value, p_value, std_err = linregress(x_values, y_values)
# Save to a CSV file
regression_results = pd.DataFrame({
    'mean': [round(avg1,2)],
    'slope': [round(slope,4)],
    'intercept': [round(intercept,2)],
    'r_value': [round(r_value,4)],
    'p_value': [round(p_value,4)],
    'std_err': [round(std_err,4)]
})
regression_results.to_csv('linear_regression_results.csv', index=False)

linear_function = lambda x: slope * x + intercept
predicted_y = linear_function(x_values)
plt.plot(x_values, predicted_y, label='Fitting Slope ~ '+'%.3f' %slope,linestyle='-', color='red')


#plt.plot(data.iloc[:, 3*(NkN+1)], 'green', label='Core ~ '+ '%.1f' % avgCore)
#plt.plot(data.iloc[:, 4*(NkN+1)], 'yellow', label='Core1 ~ '+ '%.1f' % avgCore1)
if NkN > 1:
    for i in range(NkN):
        plt.scatter(range(len(data.iloc[:, 1])),data.iloc[:, i+1], label='Pro'+str(i+1)+ ' in largest ~ '+ '%.2f' % DropProteinN[i],s=0.1)
plt.scatter(range(len(data.iloc[:, 1])),data.iloc[:, NkN+1],color='blue', label='2nd Largest ~ '+ '%.2f' % avg2,s=0.1)

plt.title(str(Nk) + ' Largest Cluster(t)') 
plt.xlabel('Simulation Progress')
plt.ylabel('Largest Cluster Size')
plt.legend(loc='lower right')
plt.savefig('Largest_Cluster.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()



#plot Normalized Largest Cluster Size
data=pd.read_csv('A_LargestCluster.txt',engine='python',sep='\s+',  header=None)
n_avg=(len(data.iloc[:,0])-1)//2 #Use last 50% as average value


t=np.linspace(0,1,len(data.iloc[:, 1]))
ClusterFrac_mRNA=data.iloc[:, 1]/Nk[0]
ClusterFrac_Pro=data.iloc[:, 2]/Nk[1]
df=pd.DataFrame({'t':t,'ClusterFrac_mRNA':ClusterFrac_mRNA,'ClusterFrac_Pro':ClusterFrac_Pro})
df.to_csv('Largest_Cluster_Normalized.csv',index=None)

plt.figure(figsize=(6.3, 3.6))
if Nk[0]!= 0:
    plt.scatter(t,100*ClusterFrac_mRNA, color='tab:red',label='Polysome: N~' + '%.2f' % (DropProteinN[0]) + ', Frac ~ ' + '%.2f' % (DropProteinN[0]/Nk[0]),s=0.1)
plt.scatter(t,100*ClusterFrac_Pro, color='tab:blue',label='Protein: N~' + '%.2f' % (DropProteinN[1]) + ', Frac ~ ' + '%.2f' % (DropProteinN[1]/Nk[1]),s=0.1)


# plt.title(f'{Nk[0]}, {Nk[1]}, Normalized Largest Cluster(t)') 
plt.xlabel('Simulation Progress')
plt.ylabel('Largest Cluster Proportion (%)')
plt.ylim(0,100)
plt.legend(loc='lower right')
plt.savefig('Largest_Cluster_Normalized.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

#plot mRNA-nucleation based Largest Cluster Size
if Nk[0]!= 0:
    plt.figure(figsize=(6.3, 3.6))
     # Split into two parts based on threshold 0.2
    mask_high = ClusterFrac_mRNA > 0.2
    mask_low = ~mask_high  # Equivalent to ClusterFrac_Pro <= 0.2

    # Plot high values (blue)
    plt.scatter(t[mask_high], 100 * ClusterFrac_Pro[mask_high], color='tab:red', s=0.1, label='Polysome  co-localized')
    
    # Plot low values (red)
    plt.scatter(t[mask_low], 100 * ClusterFrac_Pro[mask_low], color='tab:blue', s=0.1, label='Polysome not localized')


    #plt.scatter(t,100*ClusterFrac_Pro, color='tab:red',label='Protein: N~' + '%.2f' % (DropProteinN[1]) + ', Frac ~ ' + '%.2f' % (DropProteinN[1]/Nk[1]),s=0.1)


    # plt.title(f'{Nk[0]}, {Nk[1]}, Normalized Largest Cluster(t)') 
    plt.xlabel('Simulation Progress')
    plt.ylabel('Largest Cluster Proportion (%)')
    plt.legend(loc='lower right')
    plt.ylim(0,100)
    plt.savefig('Largest_Cluster_Normalized_mRNA.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
    plt.close()



# Load the data
data = pd.read_csv('A_Energy.txt', engine='python', sep='\s+', header=None)
normalized_time=np.linspace(0, 1, len(data.iloc[:, 0]))
energy=data.iloc[:, 0] / 1000
# Main figure
plt.figure(1,figsize=(6.3, 3.6))
plt.plot(normalized_time, energy, color='black', linestyle=':', label='Energy')
plt.xlim(0,1)
plt.xlabel('Relative simulation time')
plt.ylabel('Energy ($k_B$T)')


# Inset figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create an inset axis in the upper right
main_ax=plt.gca()
inset_ax1 = inset_axes(main_ax, width="45%", height="60%", loc='upper right', borderpad=1.5)
# Plot the full data in the inset
inset_ax1.plot(normalized_time, energy, color='black', linestyle=':')
# Customize the inset axis
inset_ax1.yaxis.set_label_position("left")
inset_ax1.yaxis.tick_left()
inset_ax1.set_xlim(left=0.06, right=1)  # Focus on the region of interest
#inset_ax1.set_ylim(-56000, -54000)  # Adjust limits for clarity
inset_ax1.set_ylim(min(energy[int(len(energy) * 0.1):])-500, max(energy[int(len(energy) * 0.1):])+500)  # Adjust limits for clarity
inset_ax1.set_xlabel('Relative simulation time', fontsize=8)
#inset_ax1.set_ylabel('Energy ($k_B$T)', fontsize=8)
inset_ax1.tick_params(axis='both', labelsize=8)

# Inset figure 2 (upper left)
# Create another inset axis in the upper left
plt.figure(1)
inset_ax2 = inset_axes(main_ax, width="30%", height="60%", loc='upper left', borderpad=1.5)
# Plot the first 10% of data points in the second inset
inset_ax2.plot(normalized_time, energy, color='black', linestyle=':')
# Customize the second inset axis
inset_ax2.set_xlim(0, 0.06)  # Focus on the first 10% of the time
#inset_ax2.set_ylim(energy_inset2.min() * 0.9, energy_inset2.max() * 1.1)  # Adjust limits for clarity
inset_ax2.yaxis.set_label_position("right")
inset_ax2.yaxis.tick_right()
inset_ax2.set_xlabel('Relative simulation time', fontsize=8)
inset_ax2.set_ylabel('Energy ($k_B$T)', fontsize=8)
inset_ax2.set_xticks([0, 0.05])
inset_ax2.tick_params(axis='x', labelsize=8)
inset_ax2.tick_params(axis='y', labelsize=8, direction='in', pad=-37)

#plt.legend(loc='upper right')
# Save the plot
plt.savefig('Energy.svg', format='svg', dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()



print('We are assuming sphrical droplets')
if NkN != 2:
    raise ValueError("NkN is not equal to 2, you need new analysis code, bro")

dfs = []  
for i in range(1, NkN + 1):
    filename = f'A_ClusterRadDist_Protein{i}_CenteredOnLargestCluster.txt'
    df = pd.read_csv(filename, sep='\s+', header=None)
    dfs.append(df)

n_avg=len(dfs[0])//2 #Use last 50% data to get mean
n_columns = dfs[0].shape[1]

beads_RadHist=[]
for j in range(NkN):
    beads_mean=[]
    for i in range(n_columns):
        beads_mean.append(dfs[j].iloc[-n_avg:,i].mean())
    beads_RadHist.append(beads_mean)

l=[]
hist_b1, hist_p1=[],[]
hist_b2, hist_p2=[],[]
for i in range(RD):
    l.append(i+1)
    hist_b1.append(beads_RadHist[0][i])
    hist_p1.append(hist_b1[i]/Vk[0])
    hist_b2.append(beads_RadHist[1][i])
    hist_p2.append(hist_b2[i]/Vk[1])

#Plot raw distribution
plt.plot(l,hist_b1,'r',label='mRNA Beads')
plt.plot(l,hist_b2,'b',label='Protein Beads')
plt.legend()
plt.xlabel('r')
plt.ylabel('#Beads')
plt.savefig('RadHist_Beads.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

plt.plot(l,hist_p1,'r',label='mRNA')
plt.plot(l,hist_p2,'b',label='Protein')
plt.legend()
plt.xlabel('r')
plt.ylabel('#Protiens')
plt.savefig('RadHist_Proteins.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

avg_shell_counts=[4.05, 29.44, 79.33, 154.92, 255.94, 380.12, 532.11, 708.0, 909.99, 1134.83, 1384.85, 1662.69, 1967.59, 2289.03, 2644.05, 3018.7, 3421.68, 3851.4, 4298.36, 4782.71, 5282.64, 5807.18, 6362.52, 6942.24, 7541.63, 8175.35, 8825.36, 9505.09, 10206.8, 10937.29, 11691.31, 12469.89, 13270.47, 14106.14, 14958.2, 15836.33, 16743.36, 17673.71, 18630.47, 19602.16, 20612.17, 21643.19, 22706.69, 23772.47, 24888.97, 26017.45, 27173.33, 28346.34, 29567.4, 30333.63, 30634.9, 30410.35, 29678.41, 28915.98, 28077.29, 27195.65, 26266.5, 25291.33, 24249.3, 23185.04, 22040.26, 20867.52, 19630.65, 18348.91, 17025.73, 15631.75, 14202.75, 12727.21, 11186.21, 9615.05, 8104.49, 6729.85, 5586.29, 4632.41, 3836.97, 3146.08, 2543.88, 2023.73, 1576.99, 1201.5, 872.34, 610.96, 394.05, 239.1, 117.65, 48.0, 13.19]
avg_shell_counts_2=[4.1, 29.13, 79.6, 154.52, 255.7, 380.81, 532.85, 707.18, 908.44, 1135.45, 1386.3, 1664.29, 1961.32, 2293.21, 2642.68, 3019.52, 3423.18, 3848.01, 4304.23, 4777.16, 5281.23, 5812.12, 6362.18, 6940.77, 7543.31, 8172.38, 8825.1, 9504.47, 10207.39, 10938.0, 11689.27, 12470.65, 13273.06, 14107.97, 14953.62, 15838.13, 16745.81, 17668.49, 18629.35, 19605.31, 20621.71, 21633.98, 22703.6, 23776.01, 24889.07, 26013.32, 27171.24, 28364.19, 29553.41, 30787.92, 32050.22, 33331.87, 34638.65, 35963.04, 37328.29, 38707.21, 40120.08, 41550.2, 43003.79, 44489.93, 45991.05, 47529.18, 49097.23, 50668.12, 52281.11, 53917.79, 55565.18, 57258.53, 58959.97, 60714.13, 62445.34, 64252.66, 66044.8, 67891.41, 69752.74, 71624.01, 73546.14, 75472.71, 77452.73, 79408.31, 81433.35, 83480.43, 85530.59, 87609.68, 89726.56, 91867.61, 94028.59, 96208.61, 98428.09, 100661.91, 102916.45, 105212.49, 107513.05, 109865.67, 112225.17, 114598.87, 117034.49, 119448.59, 121935.83, 123447.21, 124073.72, 123709.54, 122362.24, 120969.68, 119480.1, 118003.62, 116416.43, 114838.6, 113160.24, 111452.39, 109701.68, 107880.26, 106025.0, 104112.97, 102165.04, 100156.84, 98065.81, 95973.68, 93812.14, 91605.43, 89343.97, 87017.62, 84681.97, 82234.32, 79775.45, 77288.59, 74709.89, 72105.1, 69420.89, 66726.62, 63962.15, 61121.43, 58286.9, 55354.56, 52402.25, 49359.66, 46319.71, 43194.91, 40030.69, 36809.78, 33603.77, 30503.07, 27699.06, 25209.68, 23015.74, 21013.45, 19163.6, 17471.2, 15898.8, 14421.38, 13046.3, 11787.73, 10599.19, 9485.14, 8451.55, 7511.0, 6613.83, 5786.04, 5036.57, 4347.42, 3700.59, 3109.26, 2592.96, 2129.36, 1698.99, 1314.57, 1005.09, 732.16, 502.87, 313.43, 179.65, 85.6, 30.69, 7.27]
l=[]
rho_b1, rho_p1, rho_n1=[],[],[] #beads, protein, normalized
rho_b2, rho_p2, rho_n2=[],[],[]

if len(beads_RadHist[0]) <100:
    for i in range(RD):
        l.append(i)
        if Nk[0]!= 0:
            rho_b1.append(beads_RadHist[0][i]/avg_shell_counts[i])
            rho_p1.append(rho_b1[i]/Vk[0])
        rho_b2.append(beads_RadHist[1][i]/avg_shell_counts[i])
        rho_p2.append(rho_b2[i]/Vk[1])

else: 
    for i in range(RD):
        l.append(i)
        if Nk[0]!= 0:
            rho_b1.append(beads_RadHist[0][i]/avg_shell_counts_2[i])
            rho_p1.append(rho_b1[i]/Vk[0])
        rho_b2.append(beads_RadHist[1][i]/avg_shell_counts_2[i])
        rho_p2.append(rho_b2[i]/Vk[1])


#get normalized distribution
for i in range(RD):
    if Nk[0]!= 0:
        rho_n1.append(rho_b1[i]/sum(rho_b1))
    rho_n2.append(rho_b2[i]/sum(rho_b2))


if Nk[0]!= 0:
    plt.plot(l,rho_b1,'r',label='mRNA Beads',linewidth=3)
plt.plot(l,rho_b2,'b',label='Protein Beads',linewidth=3)
plt.legend()
plt.xlabel('r')
plt.ylabel('#Beads/#Sites')
plt.savefig('RDF_Beads.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

if Nk[0]!= 0:
    plt.plot(l,rho_p1,'r',label='mRNA',linewidth=3)
plt.plot(l,rho_p2,'b',label='Protein',linewidth=3)
plt.legend()
plt.xlabel('r')
plt.ylabel('#Protiens/#Sites')
plt.savefig('RDF_Proteins.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

plt.figure(figsize=(5, 4.5))
l_µm = [x * 0.003 for x in l]

# --- left axis ---
fig, ax1 = plt.subplots(figsize=(5, 5))

ax1.plot(l_µm, [x * C_µM/1000 for x in rho_p2], 'tab:blue', label='Protein', linewidth=3,zorder=2)
ax1.set_xlim(0, 0.26)
ax1.set_ylim(0, 2)
ax1.set_xlabel('Radial distance from condensate center (µm)', fontsize=12.5)
ax1.set_ylabel('Protein concentration (mM)', fontsize=14, color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# --- right axis ---
if Nk[0] != 0:
    ax2 = ax1.twinx()
    ax2.set_ylim(0, 0.1)
    ax2.plot(l_µm, rho_n1, 'tab:red', label='Polysome', linewidth=3,zorder=1)
    ax2.set_ylabel('Polysome density (normalized)', fontsize=14, color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
        # ensure left axis (ax1) is on top of right axis
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)  # make background transparent so both visible

# --- legends ---
lines1, labels1 = ax1.get_legend_handles_labels()
if Nk[0] != 0:
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=12)
else:
    ax1.legend(loc='upper right', fontsize=12)

plt.savefig('RDF_Normalized.svg', format='svg', dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()


# Remove the first three data points from each curve
l_trimmed = l[3:]
if Nk[0]!= 0:
    rho_p1_trimmed= rho_p1[3:]
rho_p2_trimmed = rho_p2[3:]

# Fit the trimmed data to the sigmoid function

params2, params2_covariance = curve_fit(sigmoid, l_trimmed, rho_p2_trimmed)

# Extract the fitted parameters

a2_fit, b2_fit, c2_fit = params2
if a2_fit < 1e-3:
    a2_fit=0

# Generate the fitted curve using the fitted parameters

y2_fit = sigmoid(l, a2_fit, b2_fit, c2_fit)

# Plot the original data and the fitted curve
if Nk[0]!= 0:
    plt.scatter(l, rho_p1,label='mRNA')
plt.scatter(l, rho_p2,label='Protein')


plt.plot(l, y2_fit, 'b', label='Protein_fit')

plt.legend()
plt.xlabel('r')
plt.ylabel('#Beads/#Sites')
plt.savefig('RDF_Proteins_Fitted.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.03)
plt.close()

radial_density_fitting = pd.DataFrame(columns=['Protein', 'Function', 'a', 'b', 'c'])    
radial_density_fitting.loc[0]=['Protein','a/(1+exp(b * (x - c)))','{:.2e}'.format(a2_fit), '{:.2e}'.format(b2_fit), round(c2_fit,1)]
radial_density_fitting.to_csv('RDF_Proteins_Fitting.csv', index=False)


V_c2=4*np.pi*c2_fit**3/3#size of dimer

#Create a DataFrame with the droplet info
droplet_info = pd.DataFrame(columns=['Phase', 'Size', 'Volume', 'N', 'N_mRNA', 'N_Protein', '<mRNA>', '<Protein>'])
droplet_info.loc[0] = ['Dense_a', f'{c2_fit:.1f}', f'{V_c2:.1f}', avg1,  DropProteinN[0], DropProteinN[1],
                       'N/A', f'{a2_fit:.2e}']
droplet_info.loc[1] = ['Dense_c', f'{c2_fit:.1f}', f'{V_c2:.1f}', avg1, DropProteinN[0], DropProteinN[1],
                       f'{DropProteinN[0]/V_c2:.2e}', f'{DropProteinN[1]/V_c2:.2e}']
droplet_info.loc[2] = ['Dilute_c', f'{BoxSize}', f'{V_B -V_c2:.1f}', f'{sum(Nk)-avg1:.1f}', f'{Nk[0]-DropProteinN[0]:.2f}', f'{Nk[1]-DropProteinN[1]:.2f}',
                       f'{(Nk[0]-DropProteinN[0])/(V_B-V_c2):.2e}', f'{(Nk[1]-DropProteinN[1])/(V_B-V_c2):.2e}']

#Save the output to CSV files
droplet_info.to_csv('Droplet_info.csv', index=False)

# C_µM = 10**7/(6.02*27)  # Extra 10^6 is to convert M to µM

# Create a DataFrame for the Phase Diagram
Phase_Diagram = pd.DataFrame(columns=['Type', '#mRNA', '#Protein', 'Dense%', '[Protein-]', '<mRNA>', '<Protein>', '[mRNA+]', '[Protein+]','P(Dense/Dilute)', 'N_mRNA', 'N_Protein'])
Phase_Diagram.loc[0] = ['[P-](a_fit)', Nk[0], Nk[1], round(avg1/sum(Nk), 3), f'{a2_fit:.2e}', f'{(Nk[0]/V_B):.2e}', f'{(Nk[1]/V_B):.2e}', 
                        f'{((Nk[0]-DropProteinN[0])/(V_B-V_c2)):.2e}', f'{((Nk[1]-DropProteinN[1])/(V_B-V_c2)):.2e}',f'{a2_fit/((Nk[1]-DropProteinN[1])/(V_B-V_c2))}',DropProteinN[0], DropProteinN[1],]
Phase_Diagram.loc[1] = ['a_fit /µM', Nk[0], Nk[1], round(avg1/sum(Nk), 3), f'{(C_µM*a2_fit):.2e}', f'{(C_µM*Nk[0]/V_B):.2e}', f'{(C_µM*Nk[1]/V_B):.2e}', 
                        f'{(C_µM*(Nk[0]-DropProteinN[0])/((V_B-V_c2))):.2e}', f'{(C_µM*(Nk[1]-DropProteinN[1])/(V_B-V_c2)):.2e}',f'{a2_fit/((Nk[1]-DropProteinN[1])/(V_B-V_c2))}',DropProteinN[0], DropProteinN[1],]


# Save the output to CSV files
Phase_Diagram.to_csv('Phase_Diagram.csv', index=False)

