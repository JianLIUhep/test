#!/usr/bin/env python
#
# download raw test data files or not (depending on the moon position)

from __future__ import print_function, unicode_literals
import hashlib
import io
import os
import os.path
import urllib
import sys
import tarfile
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

BASE_SOURCE = 'https://project-corryvreckan.web.cern.ch/project-corryvreckan/data/'
BASE_TARGET = 'data/'
# dataset names and corresponding checksum
DATASETS = {
    'timepix3tel_ebeam120': 'a196166ea38a14bbf00c2165a9aee37c291f1201ed39fd313cc6b3f25dfa225d',
}

def sha256(path):
    """calculate sha256 checksum for the file and return in hex"""
    # from http://stackoverflow.com/a/17782753 with fixed block size
    algo = hashlib.sha256()
    with io.open(path, 'br') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            algo.update(chunk)
    return algo.hexdigest()

def check(path, checksum):
    """returns true if the file exists and matches the checksum"""
    if not os.path.isfile(path):
        return False
    if not sha256(path) == checksum:
        return False
    print('\'%s\' checksum ok' % path)
    return True

def download(name, checksum):
    # not sure if path.join works for urls
    source = BASE_SOURCE + name + '.tar.gz'
    target = os.path.join(BASE_TARGET, name + '.tar.gz')
    if not check(target, checksum):
        print('downloading \'%s\' to \'%s\'' % (source, target))
        urlretrieve(source, target)
        if not check(target, checksum):
            sys.exit('\'%s\' checksum failed' % target)

def untar(name):
    source = os.path.join(BASE_TARGET, name + '.tar.gz')
    tar = tarfile.open(source)
    tar.extractall(BASE_TARGET)
    tar.close()
    print('\'%s\' untar ok' % source)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        datasets = DATASETS
    else:
        datasets = {_: DATASETS[_] for _ in sys.argv[1:]}
    for name, checksum in datasets.items():
        download(name=name, checksum=checksum)
        untar(name)