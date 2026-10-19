import pandas as pd

df = pd.read_parquet("data/processed/pttppg/records.parquet")
print(df.all)
# print(df.head())        # first rows
# print(df.shape)         # (rows, columns)
# print(df.dtypes)        # column types
# print(df.describe())    # basic statistics for numeric columns