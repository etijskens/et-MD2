#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module et_md2.verletlist.impl_cpp.
"""

import sys
sys.path.insert(0,'.')

import pytest


import numpy as np

from et_md2.atoms import Atoms
import et_md2.verletlist

def test_vl():
    natoms = 10
    VerletList = et_md2.verletlist.set_impl('cpp')
    vlist = VerletList(natoms,1.0);

    with pytest.raises(RuntimeError):
        vlist.has(natoms+1,9)
    assert not vlist.has(natoms-2,natoms-1)

    for i in range(0,natoms-1):
        for j in range(i+1,natoms):
            print(f'adding {(i,j)}')
            vlist.add(i,j)
            assert vlist.has(i,j)

    vlist.print()
    vlist.linearise(True)
    vlist.print()

    for i in range(0,natoms-1):
        for j in range(i+1,natoms):
            # print((i,j))
            assert vlist.has(i,j)



#===============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (normally all tests are run with pytest)
#===============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_vl

    print(f"__main__ running {the_test_you_want_to_debug} ...")
    the_test_you_want_to_debug()
    print('-*# finished #*-')
#===============================================================================
