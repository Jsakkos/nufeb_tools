def test_github_import():
    from nufeb_tools import utils
    x = utils.get_data(directory = None,test=True)

def test_local_import():
    from nufeb_tools import utils
    from pathlib import Path
    import os
    if os.path.isdir(str((Path.home()) / '.nufeb_tools' / 'data' / 'Run_26_90_83_1')):
        x = utils.get_data(directory=str((Path.home()) / '.nufeb_tools' / 'data' / 'Run_26_90_83_1'))


