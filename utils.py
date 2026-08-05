import os
import sys


def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def read_delimiter(input):
    if input:
        if 'tab' == input:
            delimiter = '\t'
        elif 'space' == input:
            delimiter = None
        elif 'comma' == input:
            delimiter = ','
        elif 'pipe' == input:
            delimiter = '|'
        else:
            delimiter = input
    else:
        delimiter = None
    return delimiter


def expand_path(path):
    """
    Create any necessary directories to ensure that the file path is valid

    :param path: a filename or directory that might or not exist
    """
    head_tail = os.path.split(path)
    if head_tail[0]:
        if not os.path.isdir(head_tail[0]):
            log('Creating directories for', head_tail[0])
            os.makedirs(head_tail[0], exist_ok=True)