import pandas as pd

df = pd.read_csv("surgeStats.csv")
df.columns = ['GD','NM','CD (J)', 'CD (G-S)']
df.index = ['Precision','Total Rides','Average Price','Iterations','Time (sec)']
df.to_latex("surge.tex")

df = pd.read_csv("surgeTaxStats.csv")
df.columns = ['CD (J)', 'CD (G-S)']
df.index = ['Precision','Total Rides','Average Price','Iterations','Time (sec)']
df.to_latex("surgeTax.tex")
