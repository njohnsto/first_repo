################
#    Validated for RS-6
#    Export DVH data to Excel 
#    This version exports data in 1 cGy dose bins, and absolute volume
#    Assumes a current plan and current beamset are loaded.
#    Allows you to select an Evaluation Sum dose if desired
#
#
#
# *****************************************************************************
# *****************************************************************************
#                  ENTER THE LOCATION OF THE RS_DVH_START.py SCRIPT HERE
SCRIPT_DIRECTORY = '//FCODC2/data/physics/RayStation/RS_DVH/'
#
import sys, os
sys.path.append(os.path.normpath(SCRIPT_DIRECTORY))
#
# ******************************************************************************
# ******************************************************************************
XLTM_DIRECTORY = SCRIPT_DIRECTORY #XLTM_DIRECTORY is the Excel template directory
DEFAULT_WBFN = XLTM_DIRECTORY + "RS_DVH.xltm"  #'Default WorkBook FileName'
#
# Store the script directory so other modules can use it
from RS_DVH_LIB import Globals

Globals.__SCRIPT_DIRECTORY = SCRIPT_DIRECTORY
Globals.__XLTM_DIRECTORY = XLTM_DIRECTORY
Globals.__DEFAULT_WBFN = DEFAULT_WBFN
#
# Keep the "start" script short, and limited to things that need no revisions.
# That way you don't have to keep uploading the script to RayStation

#Start the script
from RS_DVH_LIB import RS_DVH
