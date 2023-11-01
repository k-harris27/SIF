import json
from pathlib import Path
from typing import Dict

class ForceFields:

    # Decorator to only load FFs once.
    def _forcefield(loader):
        def inner():
            if not hasattr(loader, "ff"):
                loader.ff = loader()
            return loader.ff
        return inner

    @_forcefield
    def oplsaa() -> Dict[str,str]:
        with open(f"{ForceFields._forcefield_dir()}/oplsaa.json") as file:
            ff_text = file.read()
        return json.loads(ff_text)

    def _forcefield_dir():
        return f"{Path(__file__).resolve().parent}/type_name_lookups"
