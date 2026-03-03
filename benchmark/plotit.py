import pandas as pd
import matplotlib.pyplot as plt

def parse_time_to_minutes(time_str):
    """Converts (H:)MM:SS.ms to total minutes for plotting."""
    if pd.isna(time_str):
        return 0
    parts = str(time_str).split(':')
    if len(parts) == 3: # h:mm:ss.ms format
        mins = float(parts[0]) * 60 + float(parts[1])
        secs = float(parts[2])
        return mins + (secs / 60)
    elif len(parts) == 2: # mm:ss.ms format
        mins = float(parts[0])
        secs = float(parts[1])
        return mins + (secs / 60)
    return 0

# 1. Load Data
df = pd.read_csv('benchmark_results.csv')

# 2. Process Data
df['Temp_Ecoule_Min'] = df['Temp_Ecoule'].apply(parse_time_to_minutes)

# 3. Create Subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Définition des styles par outil
tools = df['Outil'].unique()
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
markers = ['o', 's', '^', 'D']

# Tracer les courbes pour chaque outil
for i, tool in enumerate(tools):
    tool_data = df[df['Outil'] == tool]

    # Top Plot: Elapsed Time vs. Target Size
    ax1.plot(tool_data['Taille_Cible_Go'], tool_data['Temp_Ecoule_Min'],
             marker=markers[i % len(markers)], color=colors[i % len(colors)],
             linewidth=2, label=tool)

    # Bottom Plot: Max RAM vs. Target Size
    ax2.plot(tool_data['Taille_Cible_Go'], tool_data['RAM_Max_Mo'],
             marker=markers[i % len(markers)], color=colors[i % len(colors)],
             linewidth=2, label=tool)

# Formatting Top Plot
ax1.set_ylabel('Elapsed Time (Minutes)', fontsize=11)
ax1.set_title('Benchmark Performance vs Target Input Size', fontsize=14, pad=15)
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.legend(title="Tool", title_fontsize='11', fontsize='10')

# Formatting Bottom Plot
ax2.set_xlabel('Input Size (GB)', fontsize=11)
ax2.set_ylabel('Max RAM (MB)', fontsize=11)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.legend(title="Tool", title_fontsize='11', fontsize='10')

# 4. Final Formatting
plt.tight_layout()
plt.savefig('benchmark_graphs.svg')
print("Plot successfully generated and saved as 'benchmark_graphs.svg'")