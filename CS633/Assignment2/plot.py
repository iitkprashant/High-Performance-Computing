import pandas as pd
import matplotlib.pyplot as plt
# Create a times with 4 columns using prutor assignment portal
data = {'4096_8_0': [4.7653,4.6036,4.6053,4.6022,4.7703],
        '4096_8_1': [4.5988,4.6098,4.6121,4.5984,4.6150],
        '4096_12_0': [4.6104,4.6104,4.6309,4.6428,4.6204],
        '4096_12_1': [4.7724,4.6067,4.6166,4.6034,4.6175],
        '8192_8_0': [18.3438,18.3277,18.3951,18.4000,18.3379],
        '8192_8_1': [18.2971,18.3010,18.3726,18.3432,18.3193],
        '8192_12_0': [18.6420, 18.3742, 18.3483, 18.3680,18.3450],
        '8192_12_1': [18.3814,18.3458,18.3277,18.3221,18.3410]}
df = pd.DataFrame(data)
plt.figure(figsize=(10,8))
plt.style.use('seaborn-whitegrid')
csv_file = 'time2.csv'
df.to_csv(csv_file, index=False)

# Plot boxplots for each column
df.boxplot()
plt.title('Time for Each Configuration')
plt.ylabel('Time $(sec)$')
plt.xlabel('$N, Proc, Leader$')
plt.grid(True)
plt.savefig('box_plot2.png', dpi = 400)
plt.show()
