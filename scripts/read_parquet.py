import pandas as pd

df = pd.read_parquet("data/processed/pttppg/records.parquet")
print(df.all)
