NUFEB Simulation Analysis
=========================

Get Simulation Data
-------------------
.. autoclass:: nufeb_tools.utils.get_data
    :members:
    :undoc-members:
    :member-order: bysource
    :show-inheritance:

Plotting
--------

Average Nutrient Concentration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: nufeb_tools.plot.average_nutrients

.. image:: _static/images/average_nutrients.png   
   :scale: 100%
   :align: center
   
.. code-block:: python
    
    from nufeb_tools import utils, plot
    import matplotlib.pyplot as plt
    import seaborn as sns
    x = utils.get_data(directory = None,test=True)
    f, ax = plt.subplots()
    sns.set_context('talk')
    sns.set_style('white')
    plot.average_nutrients(x.avg_con,'Sucrose',color='Green',legend=True)
