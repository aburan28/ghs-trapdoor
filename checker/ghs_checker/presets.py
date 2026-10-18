import json
from importlib.resources import files


def list_presets():
    return json.loads(files('ghs_checker').joinpath('data/presets.json').read_text())


def get_preset(name):
    for key, preset in list_presets().items():
        if name.lower() in [key.lower(), *[a.lower() for a in preset.get('aliases', [])]]:
            return preset
    raise ValueError(f"Unknown preset {name!r}; run list-presets")
