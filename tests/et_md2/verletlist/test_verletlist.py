#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""Tests for sub-module et_md2.verletlist."""

# import pytest
import sys
sys.path.insert(0,'.')


from et_md2.atoms import Atoms
import et_md2.verletlist
VerletList = et_md2.verletlist.set_impl(impl='py')
# from et_md2.grid import Grid

import numpy as np



def test_build_simple_1():
    """"""
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')

        vl = VerletList(cutoff=2.5)
        vl.build_simple(r, keep2d=True)
        print(vl)
        assert     vl.has((0, 1))
        assert     vl.has((0, 2))
        assert not vl.has((0, 3))
        assert not vl.has((0, 4))
        assert not vl.has((1, 0))
        assert     vl.has((1, 2))
        assert     vl.has((1, 3))
        assert not vl.has((1, 4))
        assert not vl.has((2, 0))
        assert not vl.has((2, 1))
        assert     vl.has((2, 3))
        assert     vl.has((2, 4))
        assert not vl.has((3, 0))
        assert not vl.has((3, 1))
        assert not vl.has((3, 2))
        assert     vl.has((3, 4))
        assert not vl.has((4, 0))
        assert not vl.has((4, 1))
        assert not vl.has((4, 2))
        assert not vl.has((4, 3))
        assert np.all(vl.vl_size == np.array([2,2,2,1,0]))
        assert np.all(vl.vl_list == np.array([1,2,2,3,3,4,4]))

def test_build_simple_2():
    """"""
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')
        vl = VerletList(cutoff=2.0)
        vl.build_simple(r, keep2d=True)
        print(vl)
        assert     vl.has((0, 1))
        assert     vl.has((0, 2))
        assert not vl.has((0, 3))
        assert not vl.has((0, 4))
        assert not vl.has((1, 0))
        assert     vl.has((1, 2))
        assert     vl.has((1, 3))
        assert not vl.has((1, 4))
        assert not vl.has((2, 0))
        assert not vl.has((2, 1))
        assert     vl.has((2, 3))
        assert     vl.has((2, 4))
        assert not vl.has((3, 0))
        assert not vl.has((3, 1))
        assert not vl.has((3, 2))
        assert     vl.has((3, 4))
        assert not vl.has((4, 0))
        assert not vl.has((4, 1))
        assert not vl.has((4, 2))
        assert not vl.has((4, 3))
        assert np.all(vl.vl_size == np.array([2,2,2,1,0]))
        assert np.all(vl.vl_list == np.array([1,2,2,3,3,4,4]))

def test_vl2pairs():
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')
        vl = VerletList(cutoff=2.0)
        vl.build_simple(r, keep2d=True)
        print(vl)
        print(et_md2.verletlist.vl2set(vl))

def test_build_simple_2b():
    """"""
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')
        vl = VerletList(cutoff=2.0)
        vl.build_simple(r, keep2d=True)
        # print(vl)
        pairs = et_md2.verletlist.vl2set(vl)
        expected = {(0,1),(0,2),(1,2),(1,3),(2,3),(2,4),(3,4)}
        assert pairs == expected

def test_neighbours():
    """"""
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')
        vl = VerletList(cutoff=2.0)
        vl.build_simple(r)
        print(vl.vl_size[0])
        assert vl.vl_size[0] == 2
        vl0 = vl.verlet_list(0)
        print(vl0)
        assert vl0[0] == 1
        assert vl0[1] == 2

def test_build_1():
    """Verify VerletList.build against VerletList.build_simple."""

    cutoff = 2.5
    x = np.array([0.0, 1, 2, 3, 4])
    n_atoms = len(x)
    r = np.empty((n_atoms,3))
    msg = ['x00','0x0','00x']
    for ir in range(3):
        r[:,:] = 0.0
        r[:,ir] = x
        print(f'*** {msg[ir]} ***')
        vl = VerletList(cutoff=cutoff)
        vl.build(r)
        print(vl)
        pairs = et_md2.verletlist.vl2set(vl)

        vlsimple = VerletList(cutoff=cutoff)
        vlsimple.build_simple(r)
        print(vlsimple)
        expected = et_md2.verletlist.vl2set(vlsimple)
        assert pairs == expected
        assert np.all(vl.vl_size == np.array([2,2,2,1,0]))
        assert np.all(vl.vl_list == np.array([1,2,2,3,3,4,4]))

def test_build_2():
    """Verify VerletList.build against VerletList.build_simple."""
    cutoff = 2.0
    atoms = Atoms()
    atoms.lattice_positions(upper_corner=(5,5,5))

    vl = VerletList(cutoff=cutoff)
    vl.build(atoms.r)
    print(vl)
    pairs = et_md2.verletlist.vl2set(vl)

    vlsimple = VerletList(cutoff=cutoff)
    vlsimple.build_simple(atoms.r)
    print(vlsimple)
    pairs_simple = et_md2.verletlist.vl2set(vlsimple)
    assert pairs == pairs_simple

# def _test_build_grid():
#     """Verify VerletList.build_grid against VerletList.build_simple."""
#     cutoff = 5.0
#     max_neighbours = 100
#     atoms = Atoms()
#     atoms.lattice_positions(upper_corner=(5,5,5))
#
#     # compute grid
#     the_grid = Grid(cell_size=cutoff, atoms=atoms, max_atoms_per_cell=100)
#     the_grid.build()
#
#     # build grid-based verlet list
#     vl = VerletList(cutoff=cutoff, max_neighbours=max_neighbours)
#     vl.build_grid(the_grid.atoms.r, the_grid)
#     print(vl)
#     pairs = et_md2.verletlist.vl2set(vl)
#
#     vlsimple = VerletList(cutoff=cutoff, max_neighbours=max_neighbours)
#     vlsimple.build_simple(atoms.r)
#     print(vlsimple)
#     expected = et_md2.verletlist.vl2set(vlsimple)
#     assert pairs == expected

# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (normally all tests are run with pytest)
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_build_1

    print(f'__main__ running {the_test_you_want_to_debug}')
    the_test_you_want_to_debug()
    print("-*# finished #*-")
# ==============================================================================
