from nufeb_tools import utils, plot
import matplotlib.pyplot as plt
x = utils.get_data(directory = None,test=True)
def avg_nute_plot():
    f, ax = plt.subplots()
    plot.average_nutrients(x.avg_con,'suc',ax=ax,color='Green',legend=True)
def single_cell_plot():
    f, ax = plt.subplots()
    x.single_cell_growth()
    plot.biomass_time(x.single_cell_biomass)
def overall_growth_plot():
    f, ax = plt.subplots()
    plot.overall_growth(x.biomass,ax=ax)
def growth_rate_plots():
    x.single_cell_growth()
    plot.growth_rate_div(x.single_cell_biomass)
    plot.growth_rate_time(x.single_cell_biomass)
    plot.growth_rate_mu(x.single_cell_biomass)
def plot_colonies():
    f, ax = plt.subplots()
    plot.plot_colony(x,time=100)
