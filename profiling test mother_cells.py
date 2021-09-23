from os import name
## After running, view with SnakeViz with snakeviz 'profiling_mother.prof'

def main():
    import cProfile
    import pstats
    from nufeb_tools import utils
    x = utils.get_data(directory=r'D:\runs\Run_50_50_100_1_2021-08-04_462847')
    with cProfile.Profile() as pr:
        
        import pandas as pd
        from scipy.spatial import KDTree
 
        
        df = x.positions.copy()
        df['mother_cell'] = -1
        IDs = df.loc[df.Timestep==0,'ID'].unique()
        ancestry = dict()
        _=[ancestry.update({ID:ID}) for ID in IDs]
        timesteps = sorted(df.Timestep.unique())
        for time in timesteps:
            for type_ in df.type.unique():
                ancestors = df[(df.type==type_) & (df.Timestep==time) & (df.ID.isin(ancestry))]
                arr1 = ancestors[['x','y','z']].to_numpy()
                tree1 = KDTree(arr1)
                motherless = df[(df.type==type_) & (df.Timestep==time) & ~(df.ID.isin(ancestry))]
                if not motherless.empty:
                    d, i = tree1.query(motherless[['x','y','z']].to_numpy(), k=1)
                    idx1 =motherless.index
                    a = ancestors.iloc[i,:].mother_cell.values
                    _ = [ancestry.update({id_:mother}) for id_,mother in zip(motherless.ID,a)]
                        
                        #df.loc[df.ID==id_,'mother_cell']=mother
        df.drop('mother_cell',inplace=True,axis=1)
        temp = pd.DataFrame.from_dict(ancestry,orient='index').reset_index()
        temp.columns=['ID','mother_cell']
        df = pd.merge(df,temp,on='ID')
        df['total_biomass'] = df.groupby(['mother_cell','Timestep']).cumsum()['biomass']

    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    #stats.print_stats()
    stats.dump_stats(filename='profiling_mother.prof')

if __name__ == '__main__':
    main()