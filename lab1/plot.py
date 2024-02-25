import matplotlib.pyplot as plt
import sys

file_paths = sys.argv[1:]

plots = []

for file_path in file_paths:
    x = []
    y = []

    with open(file_path, 'r') as file:
        for line in file:
            values = line.split()
            x.append(float(values[0]))
            y.append(float(values[1]))

    plots.append((x, y, file_path))

fig, ax = plt.subplots(figsize=(18, 5))

for plot_data in plots:
    x, y, file_path = plot_data
    ax.plot(x, y, label=file_path)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.grid(True)
legend = ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
plt.subplots_adjust(right=0.8)

legend.set_title('Legend', prop={'size': 12})

for label in legend.get_texts():
    label.set_fontsize(10)

plt.show()
