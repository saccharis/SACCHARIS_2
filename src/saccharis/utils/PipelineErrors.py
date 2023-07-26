###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import logging
import os
import sys


class PipelineException(Exception):
    def __init__(self, msg):
        self.msg = msg
        super().__init__()

    def __str__(self):
        return self.msg


class AAModelError(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class FileError(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class UserError(PipelineException):
    def __init__(self, msg, detail=None):
        super().__init__(msg)
        self.detail_msg = detail


class NewUserFile(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


def make_logger(name: str, log_dir: str, filename: str):
    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    logger = logging.getLogger(name)
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler(os.path.join(log_dir, filename))
    logger.setLevel(logging.DEBUG)
    if sys.gettrace():
        c_handler.setLevel(logging.DEBUG)
    else:
        c_handler.setLevel(logging.WARNING)
    f_handler.setLevel(logging.DEBUG)
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)
    return logger
