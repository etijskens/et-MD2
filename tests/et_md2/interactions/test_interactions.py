#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for sub-module et_md2.interactions."""
import sys
sys.path.insert(0,'.')

import pytest
import numpy as np

from et_md2.atoms import Atoms
import et_md2.interactions.lj as LJ
import et_md2.interactions
import et_md2.verletlist


def test_compute_forces_vl():
    atoms = Atoms(2,zero=True)
    r0 = LJ.R0
    atoms.r[:,:] = np.array( [ [-0.5*r0, 0.0, 0.0]
                             , [ 0.5*r0, 0.0, 0.0] ] )
    VerletList = et_md2.verletlist.implementation('py')
    vl = VerletList(cutoff=3*r0)
    vl.build_simple(atoms.r)
    assert vl.has((0,1))
    et_md2.interactions.compute_forces(atoms, vl, LJ.force_factor)
    print(atoms.a)
    assert np.all(atoms.a[0] == - atoms.a[1])
    for i in range(atoms.n):
        for d in range(3):
            assert atoms.a[i,d] == pytest.approx(0, abs=1e-15)


# def test_compute_forces_c_vl():
#     atoms = Atoms(2, zero=True)
#     r0 = LJ.R0
#     atoms.r[:,:] = np.array( [ [-0.5*r0, 0.0, 0.0]
#                              , [ 0.5*r0, 0.0, 0.0] ] )
#     vl = cpp.VList(atoms.n, 3*r0)
#     vl.add(0,1)
#     assert vl.has(0,1)
#     LJ.compute_forces(atoms, vl)
#     print(atoms.a)
#     assert np.all(atoms.a[0] == - atoms.a[1])
#     for i in range(atoms.n):
#         for d in range(3):
#             assert atoms.a[i,d] == pytest.approx(0, abs=1e-15)


def test_compute_interactions_vl():
    atoms = Atoms(2,zero=True)
    r0 = LJ.R0
    atoms.r[:,:] = np.array( [ [-0.5*r0, 0.0, 0.0]
                             , [ 0.5*r0, 0.0, 0.0] ] )
    VerletList = et_md2.verletlist.implementation('py')
    vl = VerletList(cutoff=3*r0)
    vl.build_simple(atoms.r)
    assert vl.has((0,1))
    epot = et_md2.interactions.compute_interactions(atoms, vl, LJ.potential)
    print(f"epot={epot}")
    assert epot == pytest.approx(-0.25, abs=1.e-15)


# def test_compute_interactions_c_vl():
#     atoms = Atoms(2,zero=True)
#     r0 = LJ.R0
#     atoms.r[:,:] = np.array( [ [-0.5*r0, 0.0, 0.0]
#                              , [ 0.5*r0, 0.0, 0.0] ] )
#     vl = cpp.VList(atoms.n, 3*r0)
#     vl.add(0,1)
#     assert vl.has(0,1)
# 
#     epot = LJ.compute_interactions(atoms, vl)
#     print(f"epot={epot}")
#     assert epot == pytest.approx(-0.25, abs=1.e-15)


# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (normally all tests are run with pytest)
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_compute_interactions_vl

    the_test_you_want_to_debug()
    print("-*# finished #*-")
# ==============================================================================
