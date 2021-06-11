#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
import time
import gzip
import logging
from multiprocessing.managers import BaseManager


class ProgressBar(object):
    """A simple progress bar with date stamp"""

    def __init__(self, width=50):
        """Init the ProgressBar object
        Paramters
        ---------
        width : integer, optional (default=50)
            The width of progress bar
        """
        self.last_x = -1
        self.width = width

    def update(self, x):
        """Update progress bar
        
        Paramters
        ---------
        x : integer, (0 <= x <= 100)
            the percentage of progress in [0, 100], if x equals input 
            from the last time, the progress bar will not be updated
        """
        assert 0 <= x <= 100
        if self.last_x == int(x):
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()
        if x == 100:
            sys.stderr.write('\n')

    def self_update(self):
        self.update(max(100, self.last_x + 1))


class ProgManager(BaseManager):
    pass
ProgManager.register('Prog', ProgressBar)


def get_logger(logger_name='logger', fname=None, verbosity=False):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    level = logging.DEBUG if verbosity else logging.INFO

    # LOG format
    fmt = "%(asctime)-15s [%(levelname)-5s] %(message)s"
    datefmt = "[%a %Y-%m-%d %H:%M:%S]"
    formatter = logging.Formatter(fmt, datefmt)

    # LOG file
    if fname is not None:
        fh = open(fname, 'w')
        fh.close()
        file_handler = logging.FileHandler(fname)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # LOG console
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def find_logger_basefilename(logger):
    log_file = None
    parent = logger.__dict__['parent']
    if parent.__class__.__name__ == 'RootLogger':
        for h in logger.__dict__['handlers']:
            if h.__class__.__name__ == 'FileHandler':
                log_file = h.baseFilename
    else:
        log_file = find_logger_basefilename(parent)
    return log_file


def to_str(bytes_or_str):
    """
    Return Instance of str
    """
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    """
    Return Instance of bytes
    """
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def load_fasta(fname, is_gz=False):
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        if to_str(line).startswith('>'):
            if seq_id is not None:
                sequences[seq_id] = seq
            seq_id = to_str(line).rstrip().lstrip('>')
            seq = ''
        else:
            seq += to_str(line).rstrip()
    sequences[seq_id] = seq
    f.close()
    return sequences


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)
