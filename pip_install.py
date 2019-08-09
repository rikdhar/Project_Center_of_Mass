

try:
    from pip import main as pipmain
except ImportError:
    from pip._internal import main as pipmain

pipmain(['install', 'biopandas'])

import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
ppdb = PandasPdb()

