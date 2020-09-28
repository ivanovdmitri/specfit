# Dmitri Ivanov <dmiivanov@gmail.com>
# Python interface for the cosmic ray flux fit c++ library.  All relevant
# classes and namespaces are accessible via PyROOT and can be used
# interactively in python by importing this module.
# Prefer: c++11 standard, and Python >=3.6, ROOT >=6
# Tested with g++ (GCC) 8.3.1, Python 3.6.8, ROOT 6.20/04
# Also successfully tested with g++ (GCC) 4.7.2, Python 2.7.3, ROOT 5.34/36

import sys
import re
import os

try:
    import ROOT
except ImportError:
    sys.stdout.write("ERROR: 'import ROOT' failed; try resinstalling ROOT\n")
    sys.stderr.write("and / or using correct versions of ROOT and Python\n")
    sys.exit(2)

if(not os.environ.get("SPECFIT")):
    os.environ["SPECFIT"] = os.path.dirname(os.path.abspath(__file__))
ROOT.gROOT.Macro(os.environ.get("SPECFIT")+"/load_specfit_cpplib.C")

hfiles = [f for f in os.listdir(os.environ["SPECFIT"]+"/inc") \
          if f.endswith(".h") and f != "specfit.h" and f != "specfitLinkDef.h" ]
attrlist = map(lambda s: s.strip(".h"), hfiles)

for attr in attrlist:
    try:
        exec("from " + ROOT.__name__ + " import " + attr)
    except ImportError:
        sys.stdout.write("ERROR: importing \'{:s}\' failed\n".format(x))
        sys.exit(2)
        
# clean up
del sys
del re
del os
