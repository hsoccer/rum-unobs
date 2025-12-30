import numpy as np
import pandas as pd


if __name__ == "__main__":
    # Make sure that RCM_multi_data.tex is in ./
    text = open("RCM_multi_data.tex").read()
    Raw_Data_path = "Raw Data"

    for idx_sub in range(1, 142):
        sub_data_list = text.split("\\begin{table}[H]\n\t\\centering\n\t\\begin{tabular}{ccccc}\n\t\t")[idx_sub][32:].split(" \\\\\n\t\t")[:-1]
        sub_data_df = pd.DataFrame([txt.split(" & ") for txt in sub_data_list])
        sub_data_df = sub_data_df.replace("-", np.nan)
        sub_data_df.to_csv(Raw_Data_path+"/participant_"+format(idx_sub, f"0{3}d")+".csv")
