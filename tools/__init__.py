from os.path import dirname, basename, isfile, join
import glob

# Find all .py files in current directory.
modules = glob.glob(join(dirname(__file__), "*.py"))

# Tell python all but this file are usable modules
__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
