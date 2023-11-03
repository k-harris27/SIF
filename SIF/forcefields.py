import json
from pathlib import Path
from typing import Dict

# Attributes are added dynamically when we load the json.
class ForceField(object):
    def __init__(self, atom_names:Dict[str,str]=None, topo_equivs:Dict[str,str]=None):
        self.atom_names = atom_names
        self.topo_equivs = topo_equivs


class ForceFields:

    # Decorator to only load FFs once and return as an object.
    def _forcefield(loader):
        def inner():
            if not hasattr(loader, "ff"):
                ff_dict = loader()
                ff_obj = ForceField(**ff_dict)
                loader.ff = ff_obj
            return loader.ff
        return inner
    
    # Decorator to create properties of the class itself
    class staticproperty(property):
        def __get__(self, owner_self, owner_cls):         
            return self.fget()

    @staticproperty
    @_forcefield
    def oplsaa() -> ForceField:
        with open(f"{ForceFields._forcefield_dir()}/oplsaa.json") as file:
            ff_text = file.read()
        return json.loads(ff_text)

    def _forcefield_dir():
        return f"{Path(__file__).resolve().parent}/type_name_lookups"
    