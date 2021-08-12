# -*- coding: utf-8 -*-

"""
Module et_md2.interactions
==========================

A submodule for computing the interactions (forces, interaction energy) from a Verlet list.

"""

import et_md2.verletlist
import et_md2.verletlist.c_vl

# impl = None

# def implementation(impl):
#     """Select the implementation.
#
#     :param str impl: 'py'|'cpp'
#     :raises ValueError: if impl unknown.
#     """
#     global VerletList
#
#     if impl == 'py':
#         VerletList = VL
#
#     elif impl == 'cpp':
#         import et_md2.verletlist.c_vl
#         VerletList = et_md2.verletlist.c_vl.VL
#
#     else:
#         raise ValueError(f'Unknown VerletList implementation: {impl}.')
#
#     return VerletList


def compute_forces(atoms, vl, ff=None):
    """

    :param atoms:
    :param vl:
    :param ff_functor: force factor functor
    :return:
    """
    r = atoms.r
    a = atoms.a
    if isinstance(vl, et_md2.verletlist.VL):
        for i in range(vl.natoms):
            o = vl.vl_offset[i]
            n = vl.vl_size[i]
            ri = r[i,:]
            ai = a[i,:]
            for k in range(o,o+n):
                j = vl.vl_list[k]
                rij = r[j,:] - ri
                rij2 = np.dot(rij,rij)
                rij *= ff(rij2)
                ai     += rij
                a[j,:] -= rij

        if atoms.m.shape[0] == 1:
            for i in range(atoms.n):
                a[i,:] /= atoms.m[0]

        elif atoms.m.shape[0] == atoms.n:
            for i in range(atoms.n):
                a[i,:] /= atoms.m[i]

    elif isinstance(vl, et_md2.verletlist.c_vl.VL):
        cpp.compute_forces(atoms.r, atoms.a, atoms.m, vl)


def compute_interactions(atoms, vl, potential=None):
    r = atoms.r
    epot = 0.0
    if isinstance(vl, et_md2.verletlist.VL):
        for i in range(vl.natoms):
            o = vl.vl_offset[i]
            n = vl.vl_size[i]
            for k in range(o,o+n):
                j = vl.vl_list[k]
                rij = r[j,:] - r[i,:]
                rij2 = np.dot(rij,rij)
                epot += potential(rij2)
    elif isinstance(vl, et_md2.verletlist.c_vl.VL):
        epot = cpp.compute_interactions(atoms.r, vl)

    return epot