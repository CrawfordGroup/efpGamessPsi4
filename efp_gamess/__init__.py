# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Plugin docstring.

"""
__version__ = "0.1"
__author__ = "J. Coleman Howard"

# Load C++ plugin
import os

import psi4

from .pymodule import *

plugdir = os.path.split(os.path.abspath(__file__))[0]
sofile = plugdir + "/" + os.path.split(plugdir)[1] + ".so"
psi4.core.plugin_load(sofile)
