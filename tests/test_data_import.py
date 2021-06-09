
def test_local_import():
  from nufeb_tools import utils
  x = utils.get_data(directory=r'data/Run_98_53_11_1')

def test_github_import():
    from nufeb_tools import utils
    x = utils.get_data(directory = None,test=True)
