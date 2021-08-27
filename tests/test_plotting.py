from nufeb_tools import utils, plot
import matplotlib.pyplot as plt
x = utils.get_data(directory = None,test=True)
def avg_nute_plot():
    f, ax = plt.subplots()
    plot.average_nutrients(x.avg_con,'suc',ax=ax,color='Green',legend=True)
def single_cell_plot():
    f, ax = plt.subplots()
    plot.biomass_time(x.positions)
def overall_growth_plot():
    f, ax = plt.subplots()
    plot.overall_growth(x.biomass,ax=ax)
def growth_rate_plots():
    plot.growth_rate_div(x.positions)
    plot.growth_rate_time(x.positions)
    plot.growth_rate_mu(x.positions)
def plot_colonies():
    f, ax = plt.subplots()
    plot.colony(x,time=100)
