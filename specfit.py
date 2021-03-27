#!/usr/bin/env python3

# Dmitri Ivanov <dmiivanov@gmail.com>

# Main module that performs the joint fit of the variety of cosmic ray spectrum results into
# broken power law functions and using energy scale correction functions if necessary.  
# The broken power law functions can be expanded to use arbitrary functions if necessary.

#  NOT WRAPPING THIS TO A MAIN FUNCTION TO BE ABLE TO
#  PROCEED WITH THE ANALYSIS INTERACTIVELY USING EITHER ROOT OR PYTHON
#  TERMINAL, IF NECESSARY, AND HAVE ACCESS TO ALL VARIABLES AT ALL TIME
#  FROM EITHER PYTHON OR ROOT

# Prefer: c++11 standard, and Python >=3.6, ROOT >=6
# Tested with g++ (GCC) 8.3.1, Python 3.6.8, ROOT 6.20/04
# Also successfully tested with g++ (GCC) 4.7.2, Python 2.7.3, ROOT 5.34/36
# (using sys.stdout, sys.stderr rather than print statements to be able to
# work with Python2 systems)


import sys
import os
from fnmatch import fnmatch
import argparse
import random
from datetime import datetime
from collections import defaultdict, OrderedDict
import numpy as np
# check for pandas
have_pandas=True
try:
    import pandas as pd
except ImportError:
    have_pandas = False
import ROOT;  ROOT.PyConfig.IgnoreCommandLineOptions = True
# if SPECFIT environmental variable was not set, set the environmental 
# variable for the main directory to the directory where this script was found
if(not os.environ.get("SPECFIT")):
    os.environ["SPECFIT"] = os.path.dirname(os.path.abspath(__file__))
from specfit_cpplib import TCRFluxFit, TCRFlux, TSPECFITF1, specfit_uti, specfit_canv
from flux_functions import FLUX_FUNCTIONS
from encorr_functions import CONSTANT_ENCORR_FUNCTIONS, NONLINEAR_ENCORR_FUNCTIONS

# date time token for creating file names for the plots that start with YYYYMMDD_HHMMSS
_dt_tok_ = "%dt"

# this flag allows running on data or pre-defined simulations, provided that they are available
__spectrum_results__= os.environ.get("spectrum_results")

# by default, use data
if not __spectrum_results__:
    __spectrum_results__ = "use_data"
else:
    if "simulation" in __spectrum_results__ and not "data" in __spectrum_results__:
        __spectrum_results__ = "use_simulation"
    elif "data" in __spectrum_results__ and not "simulation" in __spectrum_results__:
        __spectrum_results__ = "use_data"
    else:
        sys.stderr.write("ERROR: could not parse the meaning")
        sys.stderr.write(" of the environmental variable \'spectrum_results\'\n")
        exit(2)

# default choice of the flux fitting function
flux_function = FLUX_FUNCTIONS.get("fJ3B_18")

# Dictionary of available spectrum data mapped by the spectrum measurement name
AvailableSpectrumData = defaultdict(None)

# list of spectra that are used in the fit
spectrum_list = []

# energy scale correction function needs to be carefully implemented
# for the data depending on the experiment.  For simulations (default) 
# the implementation is simple and typically takes on one functional form
# as shown below.
def SetEncorrFunctions(constant_encorr_function,nonlinear_encorr_function):
    '''apply the energy scale correction functions to the simulated fluxes'''
    encorr_function = TSPECFITF1.Add("fENCORR","1.0",constant_encorr_function, nonlinear_encorr_function,1.0,1.0)
    for _,data in AvailableSpectrumData.items():
        data[2] = encorr_function

# Load all available flux results that are stored in the ASCII files and / or HDF5 
# files that have been found in some "path". 
def load_fluxes(path):
    '''load flux measurements from ASCII files and / or Pandas HDF5 files found in the path'''
    for root,dirs,files in os.walk(path):
        for f in files:
            fname=os.path.join(root,f)
            if any(fnmatch(fname,"*"+s) for s in [".h5",".hd5",".hdf5"]):
                if(not have_pandas):
                    sys.stderr.write("WARNING: not using {:s} because don't have pandas\n",fname)
                    continue
                hd5 = pd.HDFStore(fname) # compatible with Python2
                for key in hd5.keys():
                    try:
                        df = hd5[key]
                    except:
                        sys.stderr.write("WARNING (problem with pandas): failed to get {:s} from {:s} -- not using\n"\
                                             .format(key,fname))
                        continue
                    if hasattr(df,"iloc"):
                        cols=[np.array(df.iloc[:,k],dtype=np.float64) for k in range(4)]
                    else:
                        cols=[np.array(df[c],dtype=np.float64) \
                              for c in ["log10en","log10en_bsize","nevents","exposure"]]
                    name = key
                    if name[0].isdigit():
                        name[0]="_"
                    if len(name) > 1 and name[0] == "/":
                        name=name[1:]
                    for s in [".","/"," ","-"]:
                        name = name.replace(s,"_")
                    flux=TCRFlux(name,"Result {:s} from {:s}".format(name,os.path.basename(fname)))
                    flux.Load(len(cols[0]),*cols)
                    AvailableSpectrumData[name] = [flux,flux.GetTitle(),None]
                hd5.close()
            for suf in [".txt",".dat",".asc"]:
                if fnmatch(fname, "*"+suf):
                    name  = os.path.basename(fname).strip(suf)
                    title = "Result {:s} from {:s}".format(name,os.path.basename(fname)) 
                    encorr_function = None
                    AvailableSpectrumData[name] = [fname,title,encorr_function]
                    break
    

# if runing on simulated results
if __spectrum_results__ == "use_simulation":
    load_fluxes(os.environ["SPECFIT"]+"/sim")
    # default choice  - random sample of 8
    random.seed(datetime.now())
    spectrum_list = sorted(random.sample(list(AvailableSpectrumData), k=min(8,len(AvailableSpectrumData))))


# if running on data
if __spectrum_results__ == "use_data":
    load_fluxes(os.environ["SPECFIT"]+"/data")
    spectrum_list = sorted(list(AvailableSpectrumData))


# default choice of the energy correction functions
constant_encorr_function=CONSTANT_ENCORR_FUNCTIONS.get("fNOCONSTCORR")    # no constant correction
nonlinear_encorr_function=NONLINEAR_ENCORR_FUNCTIONS.get("fNONONLINCORR") # no nonlinear correction


# parse the command line arguments
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
spec_help = "Choose a comma-separated list of spectra to fit into one function from the following:\n"
for key,val in OrderedDict(sorted(AvailableSpectrumData.items(), key=lambda k: k[0])).items():
    spec_help += "{:25s} - {:s}\n".format("\'"+key+"\'",val[1])
spec_help += "Default: \'{:s}\'".format(",".join(spectrum_list))
if(__spectrum_results__ == "use_simulation"):
    spec_help += "(random {:d})".format(len(spectrum_list))
spec_help += "\n"
parser.add_argument("-spec", action = "store", dest = "spectrum_list", default = ",".join(spectrum_list), help = spec_help)
flux_fun_help = "Choose the flux fit function from the following:\n"
for key,val in FLUX_FUNCTIONS.items():
    # all functions have number of parameter method
    par_desc="{:d} parameters".format(val.GetNpar())
    # but only the ones that inherit from TBPLF1 have
    # a more specific method that gets the number of break points
    # in the broken power law function.  Treat that properly.
    try:
        nbreaks = getattr(val,"GetNbreaks")()
        par_desc +=",{:2d} breaks".format(nbreaks)
    except AttributeError:
        pass
    flux_fun_help += "{:10s} {:s}, {:.2f} - {:.2f}".format("\'"+key+"\'",par_desc,val.GetXmin(),val.GetXmax())
    if(key == flux_function.GetName()):
        flux_fun_help += "  --- Default"
    flux_fun_help += "\n"
parser.add_argument("-fun", action = "store", dest = "flux_function", default = flux_function.GetName(), help = flux_fun_help)
encorr_help = "Choose the energy scale correction functions, constant(linear) and nonlinear terms,\n"
encorr_help += "separating them by comma. The format is\n"
encorr_help += "\'CONSTANT_CORRECTION,NONLINEAR_CORRECTION\'\n"
encorr_help += "where the firt place is taken by one of the following costant(linear) corrections:\n"
for key,val in CONSTANT_ENCORR_FUNCTIONS.items():
    encorr_help += "{:15s} {:d} parameters, {:10s} {:.2f} - {:.2f}"\
        .format("\'"+key+"\'",val.GetNpar(),"\'"+val.GetExpFormula().Data()+"\',",val.GetXmin(),val.GetXmax())
    if(key == constant_encorr_function.GetName()):
        encorr_help +="  --- Default"
    encorr_help += "\n"
encorr_help += "and the second place is taken by one of the following nonlinear corrections:\n"
for key,val in NONLINEAR_ENCORR_FUNCTIONS.items():
    encorr_help += "{:20s} {:d} parameters, {:36s} {:.2f} - {:.2f}"\
        .format("\'"+key+"\'",val.GetNpar(),"\'"+val.GetExpFormula().Data()+"\',",val.GetXmin(),val.GetXmax())
    if(key == nonlinear_encorr_function.GetName()):
        encorr_help +="  --- Default"
    encorr_help += "\n"
encorr_help += "Default: \'%(default)s\'\n"
parser.add_argument("-encorr", action = "store", dest = "encorr_functions", \
                    default = ",".join((constant_encorr_function.GetName(),nonlinear_encorr_function.GetName())),\
                    help = encorr_help)
parser.add_argument("-log10en_min", "--log10en_min", action = "store", dest="log10en_min",default=18.0,help="minimum log10(E/eV), (Default:  %(default)s)")
parser.add_argument("-log10en_max", "--log10en_max", action = "store", dest="log10en_max",default=21.0,help="maximum log10(E/eV), (Default:  %(default)s)")
logEshld_default = 19.1
logEshld_fixed = None
parser.add_argument("-fix_shld", "-fix-shoulder-energy", "-fix-shoulder-log10en",\
                    action = "store",type=float,nargs="*",dest="logEshld_fixed",default=None,
                    help="Fix sholder energy at {:.2f} ? Providing a number will specify that energy"\
                    .format(logEshld_default))
parser.add_argument("-b", action = "store_true", dest="batch_mode",help="batch mode")
parser.add_argument("-q", action = "store_true", dest="quit_after_finishing",help="Quit after finishing")
parser.add_argument("-s", action = "store", dest="basename", default = None, \
                        help = "Save plots, provide \'some_basename_\'; use \'{:s}_\' to substitute date_time_ for some_basename_"\
                        .format("%"+_dt_tok_))
args = parser.parse_args()


if args.logEshld_fixed != None:
    if(len(args.logEshld_fixed) > 0):
        logEshld_fixed = args.logEshld_fixed[0]
    else:
        logEshld_fixed = logEshld_default
if logEshld_fixed != None:
    print(logEshld_fixed)

# Set the list of cosmic ray energy spectrum results to use in the fit
spectrum_list_given = [ name for name in args.spectrum_list.split(",") if len(name) ]

# check that the spectrum result exists and
# remove any duplicates while preserving the order of addition
spectrum_list=list()
spectrum_set=set()
for s in spectrum_list_given:
    if s not in spectrum_set:
        if s not in AvailableSpectrumData:
            sys.stderr.write("WARNING: spectrum {:s} does not exist\n".format(s))
            continue
        spectrum_list.append(s)
        spectrum_set.add(s)
    else:
        sys.stderr.write("WARNING: spectrum {:s} listed more than once; using it only once \n".format(s))
SpectrumFitData = { key: AvailableSpectrumData[key] for key in spectrum_list if key in AvailableSpectrumData }


# Set the flux fit function
flux_function = FLUX_FUNCTIONS.get(args.flux_function)
if flux_function == None:
    sys.stderr.write("ERROR: flux fit function \'{:s}\' doesn\'t exist!\n", args.flux_function)
    sys.stderr.write("Use only one of the following {:d} (cAsE sensitive): {:s}\n"\
                     .format(len(FLUX_FUNCTIONS)," ".join(["\'"+f+"\'" for f in FLUX_FUNCTIONS])))
    exit(2)

# Set the energy scale correction functions, constant(aka "linear") plus nonlinear contributions
encorr_functions = args.encorr_functions.split(",")
if len(encorr_functions) != 2:
    sys.stderr.write("ERROR: provide two comma separated correction function names,\n"+
                     "one for the constant function and another for the nonlinear correction\n");
    for line in encorr_help.split("\n"):
        if "%(default)s" in line:
            continue
        sys.stderr.write(line+"\n")
    
    exit(2)
constant_encorr_function,nonlinear_encorr_function = encorr_functions
constant_encorr_function = CONSTANT_ENCORR_FUNCTIONS.get(constant_encorr_function)
nonlinear_encorr_function = NONLINEAR_ENCORR_FUNCTIONS.get(nonlinear_encorr_function)
if constant_encorr_function == None:
    sys.stderr.write("ERROR: constant energy correction \'{:s}\' doesn\'t exist!\n".format(encorr_functions[0]))
    sys.stderr.write("Use only one of the following {:d} (cAsE sensitive): {:s}\n"\
                     .format(len(CONSTANT_ENCORR_FUNCTIONS)," ".join(["\'"+f+"\'" for f in CONSTANT_ENCORR_FUNCTIONS])))
    exit(2)
elif nonlinear_encorr_function == None:
    sys.stderr.write("ERROR: nonlinear energy correction \'{:s}\' doesn\'t exist!\n".format(encorr_functions[1]))
    sys.stderr.write("Use only one of the following {:d} (cAsE sensitive): {:s}\n"\
                     .format(len(NONLINEAR_ENCORR_FUNCTIONS)," ".join(["\'"+f+"\'" for f in NONLINEAR_ENCORR_FUNCTIONS])))
    exit(2)
else:
    SetEncorrFunctions(constant_encorr_function,nonlinear_encorr_function)

# Batch mode? (plots are not displayed but can be saved into files)
ROOT.gROOT.SetBatch(args.batch_mode)

# initialize the joint spectrum fitter with the chosen spectrum results
Fit = TCRFluxFit()

# fix shoulder in the fit?
if logEshld_fixed != None:
    if (flux_function.GetName() in ["fJ3B_18","fJ2B_19"]):
        ipar=flux_function.GetParNumber("logEshld")
        flux_function.SetParameter(ipar,logEshld_fixed)
        flux_function.SetParError(ipar,0)


Fit.SetFluxFun(flux_function)
for key,data in SpectrumFitData.items():
    name=key
    obj,title,fEnCorr=data
    obj = (obj,fEnCorr) if fEnCorr else (obj,)
    if type(obj[0]) == str:
        Fit.Add(name,title,*obj)
    else:
        Fit.Add(*obj)

# set the energy range and do the fit
Fit.SelectEnergyRange(float(args.log10en_min),float(args.log10en_max))

# do the fit and if it's successful, calculate statistical significance of the shoulder effect,
# and plot the results
globals()["__n_specfit_plots__"] = int(0)
if Fit.Fit():
    # do the statistical significance calculation for the shoulder feature at ~19.1
    # Null hypothesis means no shoulder feature.  Calculate how many events one would expect,
    # in the ansence of the feature, and then compare with the number of events observed
    # in the samples.  Use the summation of the Poisson distribution to determine the statistical significance
    # of seeing a certain number of events in the sample given the expectation in the absence of the feature.
    if flux_function.GetName() in ["fJ3B_18", "fJ2B_19"]:
        logEshld=flux_function.GetParameter("logEshld")
        logEgzk=flux_function.GetParameter("logEgzk")
        fJ_null=flux_function.MakeCopy("fJ_null",flux_function)
        fJ_null.SetLineStyle(9)
        fJ_null.SetLineColor(ROOT.kBlue)
        fJ_null.SetParameter("logEshld",21)
        fJ_null.SetParameter("logEgzk",21)
        fJ_null.SetRange(logEshld,logEgzk)
        Fit.SetNullFun(fJ_null)
        x=Fit.EvalNull()
        pch=specfit_uti.PoissonPchance(int(x.second), x.first,False)
        pch_sigma=specfit_uti.PoissonPchance(int(x.second),x.first,True)
        result="({:.2f} - {:.2f}) n_expect: {:.3f} n_observe: {:.0f} pchance = {:.3e} ({:.1f} sigma)".\
            format(logEshld,logEgzk,x.first,x.second,pch,pch_sigma)
        sys.stdout.write(result+"\n")
        sys.stdout.flush()

    # Show 3 types of plots for EACH individual spectrum measurement
    plot_type=["nevent,a,e1p", "j,a,e1p", "e3j,a,e1p"]
    # initialize correct number of canvases needed to fit all the plots
    specfit_canv.init_canvases(len(plot_type)*Fit.GetNfluxes(),min(8,Fit.GetNfluxes()))
    globals()["__n_specfit_plots__"] = specfit_canv.get_ncanvases()
    specfit_canv.set_log(0,1,0) # want log Y scale
    # go over each canvas and plot the corresponding results
    for iplot in range(len(plot_type)*Fit.GetNfluxes()):
        itype=iplot//Fit.GetNfluxes()
        iflux=iplot%Fit.GetNfluxes()
        flux=Fit.GetFlux(iflux)
        specfit_canv.cd(iplot+1)
        globals()["c{:d}".format(iplot+1)] = specfit_canv.get_canvas(iplot+1)
        flux.Draw(plot_type[itype])
    specfit_canv.Update()
    if specfit_canv.get_ncanvases():
        specfit_canv.cd(1)

    def get_fit_stats(fit):
        '''print out the fit statistics from the TFluxFit object'''
        log_likelihood_var_names = [ "log_likelihood", "log_likelihood_nonzero", "log_likelihood_restricted"]
        out = []
        for vn in log_likelihood_var_names:
            v=getattr(fit,vn)
            lgl=v.first
            ndof=int(v.second)-int(fit.nfitpar)
            lgl_p_dof = lgl/float(ndof) if ndof > 0 else lgl
            out.append("{:s} / ndof = {:.2f} / {:d} = {:.1f} Prob. = {:.1e}".\
                       format(vn, lgl, ndof, lgl_p_dof, ROOT.TMath.Prob(lgl,ndof)))
        return out
    sys.stdout.write("\n".join(get_fit_stats(Fit))+"\n")
    sys.stdout.flush()


def scan_parameter(ipar, npts = 40, par_lo = 0.0, par_up = 0.0, calc_deltas = True):
    npl = specfit_canv.get_ncanvases()
    if npl <  globals()["__n_specfit_plots__"] + 1:
        specfit_canv.init_canvases(1)
        npl = specfit_canv.get_ncanvases()
        globals()["c{:d}".format(npl)] = specfit_canv.get_canvas(npl)
    specfit_canv.cd(npl)
    g=Fit.scan_parameter(ipar,npts,par_lo,par_up,calc_deltas)
    g.Draw("a,L")
    specfit_canv.Update(npl)

def expand_datetime_tok(tok, cmd):
    '''expand string token into a full date time string'''
    newcmd = cmd
    # support escaping "%dt" with "\%dt"
    random.seed(datetime.now())
    tmp_tok="{:d}_{:d}"\
        .format(int((datetime.now()-datetime(1,1,1)).total_seconds()),random.randint(0,2**31))
    if("\\" + tok in cmd):
        newcmd = newcmd.replace("\\"+tok,tmp_tok)
    newcmd = newcmd.replace(tok,datetime.now().strftime("%Y%m%d_%H%M%S"))
    newcmd = newcmd.replace(tmp_tok,tok)
    return newcmd

if(args.basename):
    basename = expand_datetime_tok(_dt_tok_,args.basename)
    specfit_canv.save_plots(basename)
else:
    if args.quit_after_finishing:
        sys.stdout.write("\nTo save the plots, use \'-s\' option with the argument \'my_plots_\'\n")
    else:
        sys.stdout.write("\nTo save plots, call \'specfit_canv::save_plots(\"my_plots_\")\' in root mode\n")
        sys.stdout.write("or \'specfit_canv.save_plots(\"my_plots_\")\' in py mode\n")
    sys.stdout.write("You can also use the keyword \'{:s}_\' instead of \'my_plots_\' to have\n".format(_dt_tok_)) 
    sys.stdout.write("date and time as basename for your plots, using YYYYMMDD_HHMMSS_ format\n")

if(not args.quit_after_finishing):
    sys.stdout.write("\nTo switch between root and py command modes, type \'root\' or \'py\' and hit <ENTER>\n")
    sys.stdout.write("(to quit, type \'.q\' and press <ENTER>)\n\n")
else:
    sys.stdout.write("\'-q\' option used; quitting\n")

sys.stdout.write("To learn more about the options, run (in Linux terminal) \'{:s} --help | more\'\n\n".format(sys.argv[0]))

# quit if quit after finishing was requested
if(args.quit_after_finishing):
    sys.exit(0)


# Python2 vs Python3 difference handled
try:
    from __builtin__ import raw_input as the_raw_input
except ImportError:
    from builtins import input as the_raw_input


def preprocess_line(line):
    '''process the line before passing it to ROOT or Python'''
    return expand_datetime_tok(_dt_tok_,line)

def root_process_line(line):
    '''process the line in the ROOT terminal'''
    ROOT.gROOT.ProcessLine(line)

def py_process_line(line):
    '''process the line in the Python terminal'''
    exec(line,globals())
    
prmpts = {"ROOT" : {"line_text": "root", "line_count": 0, "line_processing_function" : root_process_line},
          "Python" : {"line_text": "py", "line_count": 0, "line_processing_function" : py_process_line}}
prmpt  = prmpts["ROOT"] # by default offer CERN ROOT prompt
exit_keys=[".q","quit","exit"]
exit_keys.extend([x + "()" for x in exit_keys]+[x + "();" for x in exit_keys])
exit_keys_ml=max(len(x) for x in exit_keys)
while True:
    try:
        line=the_raw_input(prmpt["line_text"]+" [{:d}] ".format(prmpt["line_count"]))
    except (KeyboardInterrupt, EOFError) as e:
        if(type(e) == KeyboardInterrupt):
            sys.stdout.write("\n"+repr(e)+" -- Exiting ...\n")
            sys.stdout.write("( while running, use commands:")
            sys.stdout.write((", ".join([ "\n\'"+x+"\'" if i%5 == 0 else "\'"+x+"\'" for i,x in enumerate(exit_keys)]))+"\n")
            sys.stdout.write("to quit. )\n")
            sys.stdout.write("( Also, while running, you can switch between ROOT and Python\n")
            sys.stdout.write("command modes by typing \'root\' or \'py\'. )\n")
            sys.stdout.write("-- Done\n")
        else:
            sys.stdout.write("\n")
        sys.exit(1)
    line = preprocess_line(line)
    if (line.lower()[:exit_keys_ml] in exit_keys):
        break
    if (line.lower() in ["root","rt",".rt",".root"]):
        prmpt = prmpts["ROOT"]
    elif (line.lower() in ["python","py", ".py",".python"]):
        prmpt = prmpts["Python"]
    else:
        try:
            prmpt["line_processing_function"](line)
        except Exception as e:
            sys.stderr.write(str(e)+"\n")
        prmpt["line_count"] += 1
