from os import name
## After running, view with SnakeViz with snakeviz 'profiling_mother.prof'

def main():
    import cProfile
    import pstats
    from nufeb_tools import utils
    import pandas as pd
    with cProfile.Profile() as pr:
        x = utils.get_data(test=True)
        df = x.positions
        df['mother_cell'] = df[df.Timestep==0].ID#.astype(pd.Int64Dtype())
        for time in df.Timestep.unique():
            for cell in df[df.Timestep==time].ID:
                if df[df.ID==cell].head(1).mother_cell.isna().values[0]:
                    cell_type = df[df.ID==cell].type.unique()[0]
                    dist = x.get_neighbor_distance(cell,time)
                    dist_type = dist[dist.type==cell_type]
                    neighborID = dist_type[dist_type.Distance == dist_type.Distance.min()].ID.values[0]
                    mother = df[df.ID==neighborID].head(1).mother_cell.values[0]
                    mask = df[df.ID==cell]
                    df.loc[mask.index,'mother_cell'] = mother
        df.mother_cell = df.mother_cell.astype('Int64')

    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    #stats.print_stats()
    stats.dump_stats(filename='profiling_mother.prof')

if __name__ == '__main__':
    main()