import sys
import time
import logging
from multiprocessing.managers import BaseManager

class ProgressBar(object):
    """
    A simple progressbar
    """
    def __init__(self, width=50):
        self.last_x = -1
        self.width = width

    def update(self, x, is_force=False):
        assert 0 <= x <= 100
        if self.last_x == int(x) and is_force is False:
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()

    def self_update(self):
        self.update(max(100, self.last_x + 1))

    def close(self):
        self.update(100)
        sys.stderr.write('\n')


class ProgManager(BaseManager):
    pass
ProgManager.register('Prog', ProgressBar)


def get_logger(logger_name='logger', fname=None, debugging=False):
    """
    generate logger
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    level = logging.DEBUG if debugging else logging.INFO

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
    """
    get filename for specific logger
    """
    log_file = None
    parent = logger.__dict__['parent']
    if parent.__class__.__name__ == 'RootLogger':
        for h in logger.__dict__['handlers']:
            if h.__class__.__name__ == 'FileHandler':
                log_file = h.baseFilename
    else:
        log_file = find_logger_basefilename(parent)
    return log_file