import json
from pathlib import Path
from typing import Dict

class ForceFields:

    # Decorator to only load FFs once and return as an object.
    def _forcefield(loader):
        def inner():
            if not hasattr(loader, "ff"):
                ff_dict = loader()
                loader.ff = ff_dict
            return loader.ff
        return inner
    
    # Decorator to create properties of the class itself
    class staticproperty(property):
        def __get__(self, owner_self, owner_cls):         
            return self.fget()

    @staticproperty
    @_forcefield
    def oplsaa() -> Dict:
        with open(f"{ForceFields._forcefield_dir()}/oplsaa.json") as file:
            ff_text = file.read()
        return json.loads(ff_text)

    def _forcefield_dir():
        return f"{Path(__file__).resolve().parent}/type_name_lookups"
    