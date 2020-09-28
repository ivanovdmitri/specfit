#!/usr/bin/env python3

from specfit_cpplib import TSPECFITF1
from collections import defaultdict

# Energy scale correction parameters: 
CONSTANT_ENCORR_FUNCTIONS=[TSPECFITF1("fNOCONSTCORR","0.0",17.0,21.0),
                           TSPECFITF1("fCONSTCORR","0.052",17.0,21.0),
                           TSPECFITF1("fCONSTCORRPAR","[0]",17.0,21.0,"S0","0.052","0.01")]
NONLINEAR_ENCORR_FUNCTIONS=[TSPECFITF1("fNONONLINCORR","0.0",17.0,21.0),
                            TSPECFITF1("fNONLINCORR0","(x>19.5)*0.08",17.0,21.0),
                            TSPECFITF1("fNONLINCORR1","(x>19.0)*(0.1*(x-19.0))",17.0,21.0),
                            TSPECFITF1("fNONLINCORRPAR0","(x>[0])*[1]",17.0,21.0,"logEs,S","19.5,0.08","0.1,0.01"),
                            TSPECFITF1("fNONLINCORRPAR1","(x>[0])*([1]*(x-[0]))",17.0,21.0,"logEs,slope","19.5,0.1","0.1,0.01")]
CONSTANT_ENCORR_FUNCTIONS = defaultdict(None,{f.GetName() : f for f in CONSTANT_ENCORR_FUNCTIONS})
NONLINEAR_ENCORR_FUNCTIONS = defaultdict(None,{f.GetName() : f for f in NONLINEAR_ENCORR_FUNCTIONS})
