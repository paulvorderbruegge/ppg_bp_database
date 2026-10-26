import pandas as pd

df = pd.read_parquet("data/structured/ppgbp/bp.parquet")
print(df.all)
