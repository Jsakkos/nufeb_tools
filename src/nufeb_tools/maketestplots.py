from nufeb_tools import utils, plot
import matplotlib.pyplot as plt
import seaborn as sns
# x = utils.get_data(r'D:\CADES Files\runs\Run_15_75_56_1')
x = utils.get_data(directory = None,test=True)
#%%
# f, axes = plt.subplots(ncols=3,nrows=2)
# for ax in axes.ravel():
#     x2.plot_overall_growth(ax)
f, ax = plt.subplots()
sns.set_context('talk')
sns.set_style('white')
plot.average_nutrients(x.avg_con,'Sucrose',ax=ax,color='Green',legend=True)
f.tight_layout()
f.savefig('average_nutrients.png')
# x.avg_con
# x2.avg_con
#%%
f, ax = plt.subplots()
sns.set_context('talk')
sns.set_style('white')
x.single_cell_growth()
plot.biomass_time(x.single_cell_biomass)
f.tight_layout()
f.savefig('biomass_vs_time.png')
#%%
sns.set_style('white')
sns.set_context('talk')
f, ax = plt.subplots()
plot.overall_growth(x.biomass,ax=ax)
f.tight_layout()
f.savefig('total_biomass_vs_time.png')
#%%
f = plot.growth_rate_div(x.single_cell_biomass)
f.savefig('growth_rate_div.png')
#%%
f = plot.growth_rate_time(x.single_cell_biomass)
f.savefig('growth_rate_time.png')