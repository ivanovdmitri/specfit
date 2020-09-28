#!/usr/bin/env python3

from specfit_cpplib import TBPLF1
from collections import defaultdict

# Choices of the flux fitting functions
FLUX_FUNCTIONS = [TBPLF1("fJ2B_18",2,"J",1e-30,18.0,21.0,  # [0] start below ankle, one break after 10 EeV (2 total)
                         "const,p1,p2,p3,logEank,logEgzk","2.0,-3.25,-2.7,-4.2,18.75,19.75",",".join(["0.1"]*6)),
                  TBPLF1("fJ3B_18",3,"J",1e-30,18.0,21.0,  # [1] start below ankle, 2 breaks after 10 EeV (3 total)
                         "const,p1,p2,p3,p4,logEank,logEshld,logEgzk","2.0,-3.25,-2.7,-3.0,-5.1,18.75,19.1,19.7",",".join(["0.1"]*8)),
                  TBPLF1("fJ1B_19",1,"J",1e-33,18.8,21.0,  # [0] start after ankle, one break after 10 EeV (1 total)
                        "const,p1,p2,logEgzk","2.0,-2.7,-4.2,19.75",",".join(["0.1"]*4)),
                  TBPLF1("fJ2B_19",2,"J",1e-33,18.8,21.0,  # [1] start after ankle, 2 breaks after 10 EeV (2 total)
                         "const,p1,p2,p3,logEshld,logEgzk","6.0,-2.8,-2.9,-5.1,19.1,19.7",",".join(["0.1"]*6))]
FLUX_FUNCTIONS = defaultdict(None, {f.GetName() : f for f in FLUX_FUNCTIONS})

