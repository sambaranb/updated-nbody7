# Match the nbody6xx pattern: re-export Nbody7 from this package's
# own interface module. Using a relative `.interface` import avoids
# the circular dependency with the framework-side stub at
# amuse.community.nbody7 (which itself dynamically loads us via
# load_code("nbody7", ...)).
from .interface import Nbody7

