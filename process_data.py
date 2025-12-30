import numpy as np
import pandas as pd


if __name__ == "__main__":
    DATA_DIR = "Raw Data/"

    dfs = []
    filenames = ["/participant_"+format(idx_sub, f"0{3}d")+".csv" for idx_sub in range(1, 142)]
    for filename in filenames:
        file = DATA_DIR + filename
        df = pd.read_csv(file, index_col=0).astype("Int64")
        dfs.append(df)

    df_collective_choice = sum(dfs) / (len(filenames) * 6)
    choice_sets = pd.DataFrame(np.where(df_collective_choice.notnull())).T.groupby(0).apply(lambda x: tuple(x.values[:, 1]))
    df_collective_choice.index = choice_sets
    df_collective_choice.to_csv("collective_choice_data.csv")