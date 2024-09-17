import pandas as pd

# Create a times with 4 columns using prutor assignment portal
data = {'512_5': [0.0697, 0.0598, 0.0715],
        '512_9': [0.0742, 0.0856, 0.0809],
        '2048_5': [0.7081, 0.6968, 0.6050],
        '2048_9': [0.9534, 0.8571, 0.9355]}
df = pd.DataFrame(data)

csv_file = 'time.csv'
df.to_csv(csv_file, index=False)

# Plot boxplots for each column
df.boxplot()
plt.title('Time for Each Configuration')
plt.ylabel('Time $(sec)$')
plt.xlabel('N, stencil')
plt.grid(True)
plt.savefig('box_plot.png', dpi = 200)
plt.show()
