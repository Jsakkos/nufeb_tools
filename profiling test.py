from os import name


def main():
    import cProfile
    import pstats
    from nufeb_tools import utils
    import pandas as pd
    with cProfile.Profile() as pr:
        x = utils.get_data(directory= r'D:\runs\Run_33_66_72_1_2021-06-24')
        fitness = list()
        timepoint = 30000
        for cell,type_ in zip(x.h5['id'][str(timepoint)].__iter__(),x.h5['type'][str(timepoint)].__iter__()):
            neigh = x.get_neighbor_distance(cell,timepoint)
            med = neigh.groupby(['type'])['Distance'].min()
            fitness.append([cell,type_,x.get_fitness(timepoint,cell),med[1],med[2]])
        df = pd.DataFrame(fitness,columns=['ID','Type','Fitness','Distance_1','Distance_2'])

    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    #stats.print_stats()
    stats.dump_stats(filename='profiling.prof')

if __name__ == '__main__':
    main()