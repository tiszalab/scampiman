#!/usr/bin/env python

import argparse
import sys, os
import subprocess
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
import time
import re
import logging
import shutil
# importing this packages functions from funcs/pctfuncs.py
# the try/except loop allows this tool to be run as `python pct.py` (.funcs), for development
# or be run as a true command line tool, `pct`, after being installed with pip (funcs)
try:
    from .funcs import pctfuncs as pctf
except:
    from funcs import pctfuncs as pctf


def pct():

    # I like to keep this print statement for debugging
    print(f"this script dir: {os.path.dirname(__file__) }")

    __version__ = "0.0.1"

    toolname = "pct"

    # this variable isn't used in this template, but it's useful for almost all bioinformatics
    def_CPUs = os.cpu_count()

    def_workdir = os.getcwd()

    parser = argparse.ArgumentParser(prog=f'{toolname}',
                                     description=f'{toolname}. the template tool. Visit '
                                    f'https://github.com/mtisza1/{toolname} for help.\n'
                                    f'Version {str(__version__)}')
    
    # here, a flag for the main command (calls str2bool in pctfuncs.py)
    parser.add_argument('-l', dest='LOUD', type=pctf.str2bool,
                          help='Make it loud?')  
      
    subparsers = parser.add_subparsers(dest='command')

    # subcommand sqawk
    sqawk_parser = subparsers.add_parser(
        'sqawk', 
        help=f'make a sqawking sound with the string.'
    )
    # arguments for the sqawk subcommand. 
    # notice 'thing' is in both subcommands as positional argument
    sqawk_parser.add_argument('thing', type=str, help='string to sqawk')

    sqawk_parser.add_argument("-c", "--fmt", dest="sqfmt",
                            type=str, choices=['lower', 'caps', 'title'],
                            default = 'caps',
                            help="format to sqawk input")
    
    # subcommand shuth

    shuth_parser = subparsers.add_parser(
        'shuth', 
        help=f'make a shuthing sound with the string.'
    )

    # arguments for the shuth subcommand. 
    shuth_parser.add_argument('thing', type=str, help='string to shuth')

    shuth_parser.add_argument("-s", "--dots", dest="shdots",
                            type=int,
                            default = 3,
                            help="trailing dots to add to input")

    # you need to run this line to parse the arguments
    args = parser.parse_args()


    #### define logger #####
    logger = logging.getLogger("pct_logger")
    logger.setLevel(logging.DEBUG)
    # stream gets printed to terminal
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    # file gets saved to a specified file
    file_handler = logging.FileHandler(os.path.join(def_workdir, 
                                                    f"pct.log"))
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    #########################

    if args.LOUD == True:
        soundr: str = 'LOUDLY: '
    else:
        soundr: str = 'quietly: '
    
    # using an if statement to use the different subcommands
    if args.command == "sqawk":
        pctf.sqawks(args.thing, args.sqfmt, soundr)

    if args.command == "shuth":
        pctf.shuths(args.thing, args.shdots, soundr)

# this is how to run the pct function by calling this script from the command line
if __name__ == "__main__":
    pct()
