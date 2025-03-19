import logging
import argparse

# since I'm logging, I need the same logger as in the main script
logger = logging.getLogger("pct_logger")

def prs(args: int) -> None:
    '''
    print the letter 's' any number of times (not used)
    '''
    for G in range(args):
        print(f's{G}')

def str2bool(v):
    '''
    converts bool-like strings to True/False boolean objects
    '''
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def sqawks(phrase: str, fmt: str, sound: str):
    '''
    some f-string tricks to sqawk in a variety of ways
    '''
    if fmt == "caps":
        logger.info(f'{sound}{phrase.upper()}!!!!')
    elif fmt == "lower":
        logger.info(f'{sound}{phrase.lower()}!!!!')
    elif fmt == "title":
        logger.info(f'{sound}{phrase.title()}!!!!')

def shuths(phrase: str, dots: int, sound: str):
    '''
    some f-string tricks to shuth with a variety of trailing dots
    '''
    logger.info(f"{sound}{phrase}{'.'*dots}")