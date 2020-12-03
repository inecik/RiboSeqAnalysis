import sys
import os
from shutil import which


def create_dir(*args):
    """
    Creates a directory if not exist.
    :param args: Strings to use in os.path.join
    :return: Path of created file
    """
    dir_path = os.path.join(*args)
    if not os.access(dir_path, os.W_OK) or not os.path.isdir(dir_path):  # Create directory if not exist
        os.mkdir(dir_path)
        print(f"{bcolors.WARNING}Directory created: {dir_path}{bcolors.ENDC}")
    return dir_path


def progressBarForTerminal (iteration, total, prefix ='Progress:', suffix ='', decimals = 1, barLength = 50):
    """
    This function should be called inside of loop, gives the loop's progress.
    :param iteration: It is integer. It is current iteration.
    :param total: It is integer. It is total iteration.
    :param prefix: It is string. It will be placed before progress bar.
    :param suffix: It is string. It will be placed after progress bar.
    :param decimals: It is integer. It is number of decimals in percent complete.
    :param barLength: It is integer. It is character length of bar.
    :return: It is void function. Nothing is returned.
    """
    filledLength = int(round(barLength * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


class bcolors:
    HEADER = '\033[95m\033[1m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def accept_directory(instruction):
    while True:
        path = input(instruction).strip()
        if os.path.isdir(path) and os.access(path, os.W_OK) and os.access(path, os.R_OK):
            return path
        else:
            print(f"It is not a valid directory: {path}")


def accept_file(instruction):
    while True:
        path = input(instruction).strip()
        if os.path.isfile(path) and os.access(path, os.R_OK):
            return path
        else:
            print(f"It is not a valid file: {path}")


def check_exist_package(pkg):
    try:
        assert which(pkg)
    except:  # AssertionError, TypeError, ModuleNotFoundError
        print(f"{bcolors.FAIL}{pkg} package should be installed.{bcolors.ENDC}")
        sys.exit(1325)


def check_exist_file(path):
    try:
        assert os.path.isfile(path) and os.access(path, os.R_OK)
    except:  # AssertionError, TypeError, ModuleNotFoundError
        print(f"{bcolors.FAIL}Path not found: {path}{bcolors.ENDC}")
        sys.exit(1325)


def scantree(mother_directory):
    """Recursively yield DirEntry objects for given directory."""
    for entry in os.scandir(mother_directory):
        if not entry.name.startswith('.'):
            if entry.is_dir(follow_symlinks=False):
                yield from scantree(entry.path)  # see below for Python 2.x
            else:
                yield entry


def get_files_metadata(mother_directory):
    scantree_iterator = scantree(mother_directory)
    output = list()
    for i in scantree_iterator:
        if not i.name.startswith('.') and i.is_file() and not i.is_symlink():
            stats = i.stat()
            output.append([i.path.split(mother_directory)[1], stats.st_size, stats.st_mtime, stats.st_ctime])
    return sorted(output)
