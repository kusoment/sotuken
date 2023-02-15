import pandas as pd
import matplotlib.pyplot as plt

csv_df = pd.read_csv('0.csv',encoding ='UTF8')
df_cut=csv_df.iloc[0:,[2,3]] #https://pythondatascience.plavox.info/pandas/%E8%A1%8C%E3%83%BB%E5%88%97%E3%81%AE%E6%8A%BD%E5%87%BA
x=df_cut[0]
y=df_cut[1]
print(df_cut[0])
plt.title("simu")
plt.grid
plt.scatter(df_cut[0],df_cut[1])
plt.show