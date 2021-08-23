===========
nufeb_tools
===========

|docs|  |pypi|  |tests|

Description
===========

Python-based tools and utilities for NUFEB simulations 

Getting Started
===============

Install via pip::

        pip install nufeb-tools

Generate NUFEB simulations from the CLI::

        nufeb-seed

Remove old runs::

        nufeb-clean

Get data from a simulation for analysis

.. code-block:: python

    from nufeb_tools import utils
    x = utils.get_data(test=True)

Plot the overall growth

.. code-block:: python

    from nufeb_tools import plot
    import matplotlib.pyplot as plt
    f, ax = plt.subplots()
    plot.overall_growth(x.biomass,ax=ax)

.. image:: /docs/_static/images/total_biomass_vs_time.png
   :align: center

Plot colonies based on initial seed cells

.. code-block:: python

    from nufeb_tools import utils, plot
    import matplotlib.pyplot as plt
    x = utils.get_data(directory= r'D:\runs\Run_21_18_56_1_2021-07-12')
    f,ax = plt.subplots()
    plot.colony(x,35000,colors,ax=ax)
    plt.show()

.. image:: /docs/_static/images/testcolony.png
   :align: center



.. |docs| image:: https://readthedocs.org/projects/nufeb-tools/badge/?version=latest
        :target: https://nufeb-tools.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. |pypi| image:: https://badge.fury.io/py/nufeb-tools.svg
        :target: https://badge.fury.io/py/nufeb-tools

.. |tests| image:: https://github.com/Jsakkos/nufeb-tools/actions/workflows/Test.yml/badge.svg
        :alt: Tox testing status

