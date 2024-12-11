#This code is an analysis scripte dedicated to Tetramer-Dimer situation.
#!/usr/bin/env python
# coding: utf-8
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
import matplotlib.colors as mcolors

protein_name=['Tetramer','Dimer']
def c(N,h):
    C_nM = 10**10/(6.02*27)  # Extra 10^9 is to convert M to nM
    V=50*50*h
    c=C_nM*N/V
    return c/1000 #micro M as unit

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

def smooth(dist,window):
    smoothed_list=[]
    for i in range(len(dist)-window):
        smoothed_list.append(sum(dist[i:window+i])/window)
    return smoothed_list
def initial_guess(dist,window):
    smoothed_dist=smooth(dist,window)
    c=(smoothed_dist[0]+smoothed_dist[-1])/2
    center=(len(dist)-window)//2
    a=smoothed_dist[center]-c
    for height in smoothed_dist:
        if height >= a/2+c:
            x1=smoothed_dist.index(height)
            break
    x0=center-x1
    return [a,1,x0,c]

BoxSize = get_variable('BoxSize')
BoxMid=(BoxSize[2]+1)/2
# Define the sigmoid function
def sigmoid(x, a, b, c):
    return a / (1 + np.exp(b * (x - c)))

def DouSig(x, a1, a2, b1, b2, c1, c2):
    return -a1 /(1 + np.exp(b1 * (x - c1)))+ a2 /(1 + np.exp(b2 * (x - c2)))

def SymSig(para, x):
    return -para[0]/(1+np.exp(para[1]*(x-(BoxMid-para[2]))))+para[0]/(1+np.exp(para[1]*(x-(BoxMid+para[2]))))+para[3]

# def dSymSig(para, x):
#     derivative=abs(para[0]*para[1]*(
#             (np.exp(para[1] * (x - (BoxMid-para[2]))) / (1 + np.exp(para[1] * (x - (BoxMid-para[2]))))**2) -
#             (np.exp(para[1] * (x - (BoxMid+para[2]))) / (1 + np.exp(para[1] * (x - (BoxMid+para[2]))))**2)
#         ))
#     return derivative

# def weight(para, x):
#     sum_weight = sum(1 / (1 + dSymSig(para, i)) for i in range(1, BoxSize[2] + 1))
#     return (1 / (1 + dSymSig(para, x))) / sum_weight


if BoxSize[0]==BoxSize[2]:
    if BoxSize[0]==100:
        RD=87
    elif BoxSize[0]==200:
        RD=174
    else:
        print('Your boxsize is neither 100 or 200, be careful about radial distribution!!')
        sys.exit()

# Read Keywords:
NkN = get_variable('NkN')
Nk = get_variable('Nk')
Vk = get_variable('Vk')
Vak = get_variable('Vak')

E_SelfLoop = get_variable('E_SelfLoop')
if E_SelfLoop!=0:
    # Read data from file
    with open("ASelfLoop.txt", "r") as file:
        data = [float(line.strip())/Nk[1] for line in file.readlines()]#assume the linker is the second protein in system
    # Plot the data
    plt.plot(range(1, len(data) + 1), data, marker='o',markersize=0.5, linestyle='-')
    plt.xlabel('Simulation Progress')
    plt.ylabel('#Self Loops')
    SelfLoopAvg=round(np.mean(data[-len(data)//2:]),3)
    plt.title(f'Number of self Loops, last 50% average ~ {SelfLoopAvg}')
    plt.savefig('SelfLoop.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
    #plt.show()
    plt.close()


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
C_nM = 10**10/(6.02*27)  # Extra 10^9 is to convert M to nM

#plot Number of Clusters of all data in subfolders
x=[]
cluster_hist=[]
#plot Smallest Cluster Size
data = pd.read_csv('AClusterHist.txt', engine='python', sep='\s+', header=None, names=['Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6', 'Column7', 'Column8', 'Column9', 'Column10', 'Column11', 'Column12', 'Column13', 'Column14', 'Column15', 'Column16', 'Column17', 'Column18', 'Column19', 'Column20', 'Column21'])
others = data.Column5 + data.Column6 + data.Column7 + data.Column8 + data.Column9 + data.Column10 + data.Column11 + data.Column12 + data.Column13 + data.Column14 + data.Column15 + data.Column16 + data.Column17 + data.Column18 + data.Column19 + data.Column20 + data.Column21

# Initialize variables
cluster_hist,N_hist = [],[]
x = []
n_avg = (len(data.Column1)-1)//10  # Use the last 10% of points
for i in range(21):
    cluster_hist.append(sum(data.iloc[-n_avg:, i]) / n_avg)
    N_hist.append((i+1)*sum(data.iloc[-n_avg:, i]) / n_avg)
    x.append(str(i+1))
x[20] = '21+'  # Rename the last bin to '21+'
N_hist[20]=sum(Nk)-sum(N_hist)

# Plot the histogram with a distinct separation for '21+'
fig, ax = plt.subplots(figsize=(6.3,3.6))

# Plot the bars for the first 20 bins

ax.bar(x[:20], cluster_hist[:20], color='tab:blue',label='Number of clusters of certain size')
ax.bar(x[:20], N_hist[:20], color='tab:cyan',label='Number of polymers in such cluster',zorder=0)

# Plot the last bar ('21+') with a hatch pattern
ax.bar(21.2, cluster_hist[20],width=1.5, color='tab:blue', hatch='//')  # Using hatch for distinction

# Add text labels for each bar
for a, b in zip(x[:20], cluster_hist[:20]):
    ax.text(a, b + 0.05, '%.1f' % b, ha='center', va='bottom', fontsize=7)

for a, b in zip(x[1:20], N_hist[1:20]):
    ax.text(a, b + 2.05, '%.1f' % b, ha='center', va='bottom', fontsize=6)


# Add text label for the '21+' bar
ax.text(21.2, cluster_hist[20] + 0.05, '%.1f' % cluster_hist[20], ha='center', va='bottom', fontsize=7)
ax.text(21.2, cluster_hist[20] + 5.05, '%.1f' % N_hist[20], ha='center', va='bottom', fontsize=6)

# Labels and title
ax.set_xlabel('Cluster Size')
ax.set_ylabel('Distribution')
ax.set_title(f'{Nk} Cluster Size Distribution')

# Save the figure
plt.legend()
plt.savefig('Cluster_Hist.svg', format='svg', dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()

plt.figure(figsize=(6.3, 3.6))
plt.plot(data.Column1,'g', label='Monomer ~ '+ '%.1f' % cluster_hist[0])
plt.plot(data.Column2, 'b', label='Dimer ~ '+ '%.1f' % cluster_hist[1])
plt.plot(data.Column3, 'y', label='Trimer ~ '+ '%.1f' % cluster_hist[2])
plt.plot(data.Column4, 'r', label='Tetramer ~ '+ '%.1f' % cluster_hist[3])
plt.plot(others ,'black', label='Others ~ '+ '%.1f' % sum(cluster_hist[4:]))
plt.title(str(Nk) + ' #Cluster(t)') 
plt.xlabel('Simulation Process')
plt.ylabel('Number of Clusters')
plt.legend(loc='upper right')
plt.savefig('Cluster_steps.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()

#plot Largest Cluster Size
data=pd.read_csv('ALargestCluster.txt',engine='python',sep='\s+',  header=None)
n_avg=(len(data.iloc[:,0])-1)//2 #Use last 50% as average value
avg1 = sum(data.iloc[-n_avg:, 0]) / n_avg #python starts with index 0!
avg2 = sum(data.iloc[-n_avg:, NkN+1]) / n_avg
#avgCore = sum(data.iloc[-n_avg:, 3*(NkN+1)]) / n_avg 
#avgCore1 = sum(data.iloc[-n_avg:, 4*(NkN+1)]) / n_avg

# Save to a CSV file
Mean_Cluster_Size = pd.DataFrame({
    'N':Nk[0]+Nk[1],
    'NT':Nk[0], 
    'ND':Nk[1],
    '1': [round(sum(data.iloc[-n_avg:, 0]) / n_avg,2)],
    '1T': [round(sum(data.iloc[-n_avg:, 1]) / n_avg,2)],
    '1D': [round(sum(data.iloc[-n_avg:, 2]) / n_avg,2)],
    '2': [round(sum(data.iloc[-n_avg:, 3]) / n_avg,2)],
    '2T': [round(sum(data.iloc[-n_avg:, 4]) / n_avg,2)],
    '2D': [round(sum(data.iloc[-n_avg:, 5]) / n_avg,2)],
    '3': [round(sum(data.iloc[-n_avg:, 6]) / n_avg,2)],
    '3T': [round(sum(data.iloc[-n_avg:, 7]) / n_avg,2)],
    '3D': [round(sum(data.iloc[-n_avg:, 8]) / n_avg,2)],
})
Mean_Cluster_Size.to_csv('Mean_Cluster_Size.csv', index=False)
DropProteinN= []
for i in range(NkN):
    DropProteinN.append(sum(data.iloc[-n_avg:, i+1]) / n_avg)
plt.figure(figsize=(6.3, 3.6))

plt.plot(data.iloc[:,0],'r', linestyle=':',label='Largest Cluster ~ '+ '%.2f' % avg1)

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
        plt.plot(data.iloc[:, i+1], linestyle=':',label='Pro'+str(i+1)+ ' in largest ~ '+ '%.2f' % DropProteinN[i])
plt.plot(data.iloc[:, NkN+1],'blue', linestyle=':',  label='2nd Largest ~ '+ '%.2f' % avg2)

plt.title(str(Nk) + ' Largest Cluster(t)') 
plt.xlabel('Simulation Progress')
plt.ylabel('Largest Cluster Size')
plt.legend(loc='lower right')
plt.savefig('Largest_Cluster.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()

#plot Normalized Largest Cluster Size
data=pd.read_csv('ALargestCluster.txt',engine='python',sep='\s+',  header=None)
n_avg=(len(data.iloc[:,0])-1)//2 #Use last 50% as average value
plt.figure(figsize=(6.3, 3.6))
plt.plot(np.linspace(0,1,len(data.iloc[:, 1])),data.iloc[:, 1]/Nk[0], color='tab:blue',linestyle=':',label='Tetramer: ' + '%.2f' % (DropProteinN[0]*Vak[0]))
plt.plot(np.linspace(0,1,len(data.iloc[:, 1])),data.iloc[:, 2]/Nk[1], color='tab:red',linestyle=':',label='Dimmer: ' + '%.2f' % (DropProteinN[1]*Vak[1]))

plt.xlim(left=0)
plt.ylim(bottom=0)
#plt.title(f'{Nk[0]*Vak[0]}, {Nk[1]*Vak[1]}, Largest Cluster(t)') 
plt.xlabel('Relative simulation time')
plt.ylabel('Largest Cluster Fraction')
plt.legend(loc='lower right')
plt.savefig('Largest_Cluster_Normalized.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()

# Load the data
data = pd.read_csv('AEnergy.txt', engine='python', sep='\s+', header=None)
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

#plot all molecule distribution as 10 frames

colors=['tab:blue','tab:red']
# beads_mean = np.zeros((NkN,10, BoxSize[2]))  # 10 frames, Boxsize (400) wide
# for i in range(NkN):
#     filename = f'AXYZDist_All_Protein{i+1}_Z.txt'
#     df = pd.read_csv(filename, sep='\s+', header=None)
#     n_frame = len(df) // 10  # Use 10% data as one frame
#     for i1 in range(10):
#         for j in range(BoxSize[2]):
#             beads_mean[i][i1][j] = df.iloc[n_frame*i1:n_frame*(i1+1), j].mean()

# # Convert list of numpy arrays to a single numpy array
# #beads_Dist = np.array(beads_Dist)

# l=np.linspace(1,BoxSize[2]+1,BoxSize[2])

# for i in range(NkN):
#     for j in range(10):
#         plt.plot(l,beads_mean[i][j],label=f'{protein_name[i]} {j}',color=colors[i],alpha=(j+1)/10)
# plt.legend()
# plt.xlabel('r')
# plt.ylabel('#Beads')
# #plt.title(str(Nk)+' Dist_All_Z')
# plt.savefig('Dist_All_Z_10.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
# plt.show()
# plt.close()


cmap_names = ['Reds', 'YlGn']

# Generate some mock data for the beads_mean array
beads_mean = np.zeros((NkN, 10, BoxSize[2]))  # 10 frames, Boxsize (400) wide
for i in range(NkN):
    filename = f'AXYZDist_All_Protein{i+1}_Z.txt'
    df = pd.read_csv(filename, sep='\s+', header=None)
    n_frame = len(df) // 10  # Use 10% data as one frame
    for i1 in range(10):
        for j in range(BoxSize[2]):
            beads_mean[i][i1][j] = df.iloc[n_frame*i1:n_frame*(i1+1), j].mean()

# Create the x-axis values
l = np.linspace(1, BoxSize[2]+1, BoxSize[2])

# Create a figure and subplots
fig, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

# Plotting
for i in range(NkN):
    cmap = cm.get_cmap(cmap_names[i])
    for j in range(10):
        alpha = (j + 1) / 10
        color = cmap(alpha)
        axs[i].plot(l, beads_mean[i][j], color=color)
    axs[i].set_ylabel('#Beads')
    axs[i].set_xlim((0,BoxSize[2]))
    axs[i].set_ylim(bottom=0)
    axs[i].legend([protein_name[i]])

# Create color maps for transparency levels
cmap_tetramer = cm.Reds
cmap_dimer = cm.YlGn
norm = mcolors.Normalize(vmin=1, vmax=10)
sm_tetramer = cm.ScalarMappable(cmap=cmap_tetramer, norm=norm)
sm_dimer = cm.ScalarMappable(cmap=cmap_dimer, norm=norm)
sm_tetramer.set_array([])
sm_dimer.set_array([])

# Add color bars
cbar_tetramer = plt.colorbar(sm_tetramer, ax=axs[0], orientation='vertical', pad=0.01)
cbar_tetramer.set_label('Tetramer')

cbar_dimer = plt.colorbar(sm_dimer, ax=axs[1], orientation='vertical', pad=0.01)
cbar_dimer.set_label('Dimer')

# Set labels and title
#axs[1].set_xlabel('r')
# Save and show the plot
plt.savefig('Dist_All_Z_10_frames.svg', format='svg', dpi=150, bbox_inches='tight', pad_inches=0.01)
#plt.show()
plt.close()

#plot all molecule distribution and fit
beads_Dist,protein_Dist = [],[]
for i in range(1, NkN + 1):
    filename = f'AXYZDist_All_Protein{i}_Z.txt'
    df = pd.read_csv(filename, sep='\s+', header=None)
    n_avg = len(df) // 2  # Use last 50% data to get mean
    n_columns = df.shape[1]
    beads_mean = np.zeros(n_columns)  # Initialize as a numpy array
    for j in range(n_columns):
        beads_mean[j] = df.iloc[-n_avg:, j].mean()
    beads_Dist.append(beads_mean)# for different proteins
    protein_Dist.append(beads_mean / Vk[i-1])

l=np.linspace(1,BoxSize[2]+1,BoxSize[2])

for i in range(NkN):
    plt.plot(l,beads_Dist[i],label=f'Protein{i+1} Beads')
plt.legend()
plt.xlabel('r')
plt.ylabel('#Beads')
plt.title(str(Nk)+' Dist_All_Z')
plt.savefig('Dist_All_Z.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
plt.close()

SlabSizes=[]
cDense, cDilute=[],[]
P=[]
fitting_results = pd.DataFrame(columns=['a','b','x1','x2','c'])
for i in range(NkN):  
    def error(para):
        return SymSig(para,l) - protein_Dist[i]

    para0=initial_guess(protein_Dist[i],10)
    lower_bounds = [0, 0, 0, 0]
    upper_bounds = [np.inf, np.inf, BoxMid, np.inf]
    fitting = least_squares(error, para0, bounds=(lower_bounds, upper_bounds))
    with open(f'fitting_MSD{i}.txt', 'w') as file:
        file.write(str(fitting))
    plt.figure(1)
    plt.plot((l),[SymSig(fitting.x,x) for x in l],
    color=colors[i],label=f'Fit{i+1}: H:{round(SymSig(fitting.x,BoxMid),3)} L:{round(SymSig(fitting.x,BoxSize[2]),3)}')
    P.append(SymSig(fitting.x,BoxMid)/SymSig(fitting.x,BoxSize[2]))
    cDense.append(C_nM *SymSig(fitting.x,BoxMid)/(50*50))# in Lattice unit
    cDilute.append(C_nM *SymSig(fitting.x,BoxSize[2])/(50*50))
    plt.scatter(l,protein_Dist[i], color=colors[i],s=1,label=f'Protien{i+1}_simulation')
    plt.figure(2,figsize=(6.3,3.6))
    plt.plot(l*0.003,[SymSig(fitting.x,x)*C_nM*Vak[i]/(BoxSize[0]*BoxSize[1]*1000) for x in l],
        color=colors[i],label=f'{protein_name[i]} fit')
    plt.scatter([x*0.003 for x in l],[y*C_nM*Vak[i]/(BoxSize[0]*BoxSize[1]*1000) for y in protein_Dist[i]], color=colors[i],s=1,label=f'{protein_name[i]}')
    SlabSizes.append(fitting.x[2]*2)
plt.figure(1)
plt.legend()
plt.xlim((0,BoxSize[2]))
plt.ylim(bottom=0)
plt.xlabel('r')
plt.ylabel('#Proteins')
plt.title(str(Nk)+' Dist_All_Protein_Z')
plt.savefig('Dist_All_Z_Fitted_MSD.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
#plt.show()
plt.close(1)


print(f'C dilute is {cDilute[0]},{cDilute[1]}')
if np.mean(P)>3:#this means it's phase separating
    plt.figure(2)
    print(f'this is a phase separating system as cDense/cDilute is {P[0]},{P[1]}')
    for i in range(NkN):
        #method1, direct average of dilute region.
        start=int(np.floor(BoxMid-(fitting.x[2]+34)))
        end=int(np.ceil(BoxMid+fitting.x[2]+34))
        dilute_list=np.concatenate((protein_Dist[i][:start],protein_Dist[i][end:]))
        dilute_mean=np.mean(dilute_list)
        print(f'removed region is  [{start}, {end}]')
        print(f'C dilute average is {dilute_mean*C_nM/(BoxSize[0]*BoxSize[1])}')
        cDilute[i]=dilute_mean*C_nM/(BoxSize[0]*BoxSize[1])#in nM and polymer concentration
        P[i]=cDense[i]/cDilute[i]
    plt.xlim((0,BoxSize[2]*0.003))
    plt.ylim(bottom=0)
    plt.legend(loc='upper right')
    plt.xlabel('z (μm)', fontsize=15)
    plt.ylabel('c (μM)', fontsize=15)
    # Create the inset plot
    ax_inset = inset_axes(plt.gca(), width="25%", height="65%", loc='upper left', borderpad=0.7)
    for i in range(NkN):
        ax_inset.plot([0,start*0.003],[cDilute[i]*Vak[i]/1000]*2, color=colors[i],linestyle='-')
        ax_inset.scatter([x*0.003 for x in l],[y*C_nM*Vak[i]/(BoxSize[0]*BoxSize[1]*1000) for y in protein_Dist[i]], color=colors[i],s=1,label=f'{protein_name[i]}')
    ax_inset.yaxis.set_label_position("right")
    ax_inset.yaxis.tick_right()
    ax_inset.tick_params(axis='x', labelsize=8)  # Set x-axis tick label size
    ax_inset.tick_params(axis='y', labelsize=8)  # Set y-axis tick label size
    ax_inset.set_xlim((0,start*0.003))
    ax_inset.set_ylim(0,(cDilute[0]*4+cDilute[1]*2)/1000+8)
    plt.savefig(f'Concentration_Profile_P={np.mean(P):.2f}.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
    #plt.show()
    plt.close()
else:
    plt.figure(2)
    plt.xlim((0,max(l)*0.003))
    plt.ylim(bottom=0)
    plt.legend(loc='lower right')
    plt.xlabel('z (μm)', fontsize=15)
    plt.ylabel('c (μM)', fontsize=15)
    plt.savefig(f'Concentration_Profile_P={np.mean(P):.2f}.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
    plt.close()
        #method 2, relative fitting, which is bad..
    #     def relative_error(para):
    #         return 1 - protein_Dist[i]/SymSig(para,l)
    #     #para0=initial_guess(protein_Dist[i],10)
    #     #lower_bounds = [0, 0, 0, 0]
    #     #upper_bounds = [np.inf, np.inf, BoxMid, np.inf]
    #     fitting = least_squares(relative_error, para0, bounds=(lower_bounds, upper_bounds))
    #     with open(f'fitting_rMSD{i}.txt', 'w') as file:
    #         file.write(str(fitting))
    #     plt.figure(1)
    #     plt.plot(l,[SymSig(fitting.x,x) for x in l],
    #     color=colors[i],label=f'Fit{i+1}: H:{round(SymSig(fitting.x,BoxMid),3)} L:{round(SymSig(fitting.x,BoxSize[2]),3)}')
    #     cDilute[i]=SymSig(fitting.x,BoxSize[2])/(50*50) # in Lattice unit
    #     plt.scatter(l,protein_Dist[i], color=colors[i],s=1,label=f'Protien{i+1}_simulation')

    #     plt.figure(2,figsize=(6.3,3.6))
    #     plt.plot(l*0.003,[SymSig(fitting.x,x)*C_nM*Vak[i]/(BoxSize[0]*BoxSize[1]*1000) for x in l],
    #         color=colors[i],label=f'{protein_name[i]} fit')
    #     plt.scatter(l*0.003,[y*C_nM*Vak[i]/(BoxSize[0]*BoxSize[1]*1000) for y in protein_Dist[i]], color=colors[i],s=1,label=f'{protein_name[i]}')

    # plt.figure(1)
    # plt.legend()
    # plt.xlim((0,BoxSize[2]))
    # plt.ylim(bottom=0)
    # plt.xlabel('r')
    # plt.ylabel('#Proteins')
    # plt.title(str(Nk)+' Dist_All_Protein_Z')
    # plt.savefig('Dist_All_Z_Fitted_rMSD.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
    # #plt.show()
    # plt.close(1)

    # plt.figure(2)
    # plt.xlim((0,max(l)*0.003))
    # plt.ylim(bottom=0)
    # plt.legend()
    # plt.xlabel('z (μm)')
    # plt.ylabel('c (μM)')
    # #plt.title(str(Nk)+' Concentration Profile')
    # plt.savefig('Concentration_Profile_rMSD.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
    # plt.show()
    # plt.close()
print(f'C dilute updated is {cDilute[0]},{cDilute[1]}')
print(f'Average slab height is {np.mean(SlabSizes):.2f}')
SlabSize=np.mean(SlabSizes)
V_c1=BoxSize[0]*BoxSize[1]*SlabSizes[0]
V_c2=BoxSize[0]*BoxSize[1]*SlabSizes[1]

Phase_Diagram = pd.DataFrame(columns=['Type', '#Tetramer', '#Dimer', 'Dense%', '[Tetramer-]', '[Dimer-]', '<Tetramer>', '<Dimer>', '[Tetramer+]', '[Dimer+]','[SelfLoop%]','P(Dense/Dilute)'])
Phase_Diagram.loc[0] = ['Height/nM', Nk[0], Nk[1], round(avg1/sum(Nk), 3), 
            f'{cDense[0]}', f'{cDense[1]}',
            f'{C_nM * Nk[0]/V_B}', f'{C_nM * Nk[1]/V_B}', 
            f'{cDilute[0]}', f'{cDilute[1]}', 
            SelfLoopAvg, round(np.mean(P),3)]

cT_Dense = DropProteinN[0] / V_c1
cD_Dense = DropProteinN[1] / V_c2
cT_Dilute = (Nk[0] - DropProteinN[0]) / (V_B - V_c1)
cD_Dilute = (Nk[1] - DropProteinN[1]) / (V_B - V_c2)

Phase_Diagram.loc[1] = ['Height*Vak/μM', Nk[0], Nk[1], round(avg1/sum(Nk), 3), 
            f'{cDense[0]*4/1000}', f'{cDense[1]*2/1000}',
            f'{C_nM * Nk[0]*4/(1000*V_B)}', f'{C_nM * Nk[1]*2/(1000*V_B)}', 
            f'{cDilute[0]*4/1000}', f'{cDilute[1]*2/1000}', 
            SelfLoopAvg, round(np.mean(P),3)]

Phase_Diagram.loc[2] = ['width/nM', Nk[0], Nk[1], round(avg1/sum(Nk), 3), 
                        f'{C_nM * cT_Dense}', f'{C_nM * cD_Dense}',
                        f'{C_nM * Nk[0]/V_B}', f'{C_nM * Nk[1]/V_B}', 
                        f'{C_nM * cT_Dilute}', f'{C_nM * cD_Dilute}', 
                        SelfLoopAvg, round((cT_Dense/cT_Dilute + cD_Dense/cD_Dilute) / 2,3)]

#plot cluster distribution
beads_SlabHist=[]
for i in range(1, NkN + 1):
    filename = f'AXYZDist_Cluster_Protein{i}_Z.txt'
    df = pd.read_csv(filename, sep='\s+', header=None)
    #dfs.append(df)
    n_avg=len(df)//2 #Use last 50% data to get mean
    n_columns = df.shape[1]
    beads_mean=[] #intermediate to store distribution of one bead
    for j in range(n_columns):
        beads_mean.append(df.iloc[-n_avg:,j].mean())
    beads_SlabHist.append(beads_mean)

l=range(1,BoxSize[2]+1)
for i in range(NkN):
    plt.plot(l,beads_SlabHist[i],label=f'Protein{i+1} Beads')
plt.legend()
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel('r')
plt.ylabel('#Beads')
plt.title(str(Nk)+' Dist_Cluster_Z')
plt.savefig('Dist_Cluster_Z.svg',format='svg',dpi=150, bbox_inches='tight', pad_inches=0.01)
#plt.show()
plt.close()

# Save the output to CSV files
Phase_Diagram.to_csv('Phase_Diagram.csv', index=False)