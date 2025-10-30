from pydantic import BaseModel, ConfigDict, ValidationError
from typing import Callable, Dict, Optional, Literal

# Enums for strict validation
Nature = Literal["bidirectional", "unidirectional"]
RateLaw = Literal["first", "zero"]
Regime = Literal["constant", "custom"]


class FluxSpec(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    source: str
    dest: str
    rate_constant: float
    rate_law: RateLaw = "first"
    nature: Nature = "bidirectional"


class ClearanceSpec(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    source: str
    rate_constant: float
    rate_law: RateLaw = 'first'


class DosageSpec(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    dest: str
    regime: Regime = 'constant'
    rate_constant: Optional[float] = 0.0
    dosage_func: Callable = None


class ModelSpec(BaseModel):
    compartments: Dict[str, float]
    fluxes: Optional[Dict[str, FluxSpec]] = None
    clearances: Optional[Dict[str, ClearanceSpec]] = None
    dosages: Optional[Dict[str, DosageSpec]] = None


if __name__ == "__main__":
    config = {
        "compartments": {
            "central": 22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source": "central",
                "dest": "peripheral",
                "rate_constant": 5.0,
                "nature": "bidirectional",
                "rate_law": "first"
            }
        },

        "clearances": None,

        "dosages": {
            "central_dosage": {
                "dest": "central",
                "regime": "constant",
                "rate_constant": 1.0,
            }
        }
    }

    try:
        spec = ModelSpec.model_validate(config)
    except ValidationError as e:
        print("Config invalid")
        print(e.json())
        raise
    else:
        print("Config valid. Parse:", spec)
