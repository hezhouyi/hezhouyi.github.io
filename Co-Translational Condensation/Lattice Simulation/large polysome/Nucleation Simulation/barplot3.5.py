import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Load your data
df = pd.read_csv("t_nuc_summary_SML.csv")
seq_names = ['Homo', 'M-Block','N-term', 'C-term', 'Cross']

# Define mapping of PolysomeType to label and color
label_map = {0: 'No polysome', 1: 'Small polysome', 2: 'Medium polysome', 3: 'Large polysome'}
polysome_short = {0: 'N', 1: 'S', 2: 'M', 3: 'L'}
seq_colors = ['tab:cyan','tab:red','tab:brown','tab:orange', 'tab:green']

# OPTION 1: Grey mixing approach (recommended)
def mix_with_grey(color, grey_ratio, polysome_type):
    """Mix base color with grey, more grey for lower polysome types"""
    base_rgb = mcolors.to_rgb(color)
    grey_rgb = (0.7, 0.7, 0.7)  # Light grey
    
    # Mix ratios: N=70% grey, S=50% grey, M=25% grey, L=0% grey
    mix_ratios = {0: 1, 1: 0.7, 2: 0.3, 3: 0.0}
    ratio = mix_ratios[polysome_type]
    
    mixed_rgb = tuple(base * (1 - ratio) + grey * ratio for base, grey in zip(base_rgb, grey_rgb))
    return mixed_rgb

# OPTION 2: Brightness adjustment (alternative)
def adjust_brightness(color, polysome_type):
    """Adjust brightness: darker for higher polysome types"""
    base_rgb = mcolors.to_rgb(color)
    # Brightness factors: N=40%, S=60%, M=80%, L=100%
    brightness_factors = {0: 0.4, 1: 0.6, 2: 0.8, 3: 1.0}
    factor = brightness_factors[polysome_type]
    
    adjusted_rgb = tuple(c * factor for c in base_rgb)
    return adjusted_rgb

# OPTION 3: Hatching patterns (alternative)
hatch_patterns = {0: '...', 1: '///', 2: 'xxx', 3: None}  # None = solid

# Sort and prepare data
df['Label'] = df['PolysomeType'].map(label_map)

# Set up bar positions
n_seq = df['Seq'].nunique()
n_type = 4
bar_height = 0.2
y = np.arange(n_seq)

fig, ax = plt.subplots(figsize=(4.54,7.8))

# Collect all y positions and labels for ticks
all_positions = []
all_labels = []

# Loop over polysome types
for i, nm in enumerate([0, 1, 2, 3]):
    sub = df[df['PolysomeType'] == nm]
    positions = y + (i - 1.5) * bar_height
    
    # Create colors - CHOOSE ONE OF THE OPTIONS BELOW:
    
    # OPTION 1: Grey mixing (most subtle and professional)
    colors = [mix_with_grey(seq_colors[seq_idx-1], 0.5, nm) for seq_idx in sub['Seq']]
    
    # OPTION 2: Brightness adjustment (uncomment to use instead)
    # colors = [adjust_brightness(seq_colors[seq_idx-1], nm) for seq_idx in sub['Seq']]
    
    # OPTION 3: Hatching (uncomment to use instead, also uncomment hatch parameter below)
    # colors = [seq_colors[seq_idx-1] for seq_idx in sub['Seq']]
    
    ax.barh(positions, sub['Mean_t_nuc'],
            xerr=sub['SE_t_nuc'],
            height=bar_height,
            color=colors,
            # hatch=hatch_patterns[nm],  # Uncomment for hatching option
            edgecolor='black',  # Add black edges for better definition
            linewidth=0.5,
            label=label_map[nm], 
            capsize=4)
    
    # Store positions and labels for y-ticks
    for j, pos in enumerate(positions):
        label = f"{polysome_short[nm]}"
        all_positions.append(pos)
        all_labels.append(label)

# Formatting
ax.set_yticks(all_positions)
ax.set_yticklabels(all_labels, fontsize=11)
ax.set_xlabel("Nucleation time $t_{nuc}$ (normalized)", fontsize=14)
ax.set_xlim(0, 1)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig("t_nuc_barplot.svg", format='svg', bbox_inches='tight', dpi=150)
plt.show()