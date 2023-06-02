import pandas as pd

#TODO Do the same for value_prop
#TODO Maye reset time step so that there is not these weird steps in the final fig
def main():
    d = pd.read_parquet("sourcesink2.parquet")
    d = d.assign(diff_L = lambda x: x.groupby(["row_id", "L"])['value'].diff().fillna(0.001)) 
    sum_diff_timestep = d.groupby(["row_id", "timestep"])['diff_L'].sum().reset_index()
    sum_diff_timestep = sum_diff_timestep[(sum_diff_timestep.diff_L >= -0.001) & (sum_diff_timestep.diff_L <= 0.0001)][['row_id', 'timestep']]
    # sum_diff_timestep = sum_diff_timestep[sum_diff_timestep.diff_L == 0 ][['row_id', 'timestep']]
    outer_join = d.merge(sum_diff_timestep, how='outer', on = ['row_id', 'timestep'], indicator=True)
    anti_join = outer_join[outer_join['_merge'] == 'left_only'].drop(columns=['_merge', 'diff_L'])

    anti_join.to_parquet("sourcesink2_simple.parquet", index=False)