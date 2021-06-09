def test_seed():
    from nufeb_tools.generate_atom import main, clean
    import sys

    main(sys.argv[1:])
    clean()

def test_seed_dev():
    import sys
    from nufeb_tools.generate_atom_dev import main, clean
    main(sys.argv[1:])
    clean()