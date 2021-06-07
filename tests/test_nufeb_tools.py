import pytest

from nufeb_tools import datafed, generate_atom, utils, plot

__author__ = "Jonathan Sakkos"
__copyright__ = "Jonathan Sakkos"
__license__ = "MIT"


#def test_get_data():
  #  """API Tests"""
  #  x = utils.get_data(directory=r'../data/Run_98_53_11_1')
  #  f, axes = plt.subplots(ncols=3,nrows=2)
   # for ax in axes.ravel():
   #     x.plot_overall_growth(ax)
  #  f, ax = plt.subplots()
  #  sns.set_context('talk')
  #  sns.set_style('white')
  #  x.plot_average_nutrients('Sucrose',color='Green')
   # assert fib(1) == 1

    #with pytest.raises(AssertionError):
        #fib(-10)


#def test_main(capsys):
  #  """CLI Tests"""
    # capsys is a pytest fixture that allows asserts agains stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
   # main(["7"])
   # captured = capsys.readouterr()
   # assert "The 7-th Fibonacci number is 13" in captured.out
