# -*- coding: utf-8 -*-

"""
Module et_md2.verletlist
======================++

A submodule for representing Verlet lists.
There is a Python and a C++ implementation.

"""
import numpy as np

VerletList = None

def set_impl(impl):
    """Select the implementation

    :param str impl: 'py'|'cpp'
    :raises ValueError: if impl unknown.
    """
    global VerletList

    if impl == 'py':
        VerletList = VL

    elif impl == 'cpp':
        import et_md2.verletlist.__impl_cpp
        VerletList = et_md2.verletlist.impl_cpp.VL

    else:
        raise ValueError(f'Unknown VerletList implementation: {impl}.')

    return VerletList


# This is the Python implementation
class VL:
    nneighbours = 20  # initial size of the Verlet lists. They grow dynamically as needed.

    def __init__(self, cutoff=1.0):
        """Verlet lists of a number of atoms.

        :param float cutoff: cutoff distance
        :param integer max_neighours: maximum size of an atom's Verlet list.
            If this number is too small build|build_simple|build_grid will
            raise IndexError.

        The initial data structure is a 2D integer numpy array. There is one
        row for each atom. Each row starts with the number of neighbours,
        followed by the atom indices of the neighbours. Thus:

        * vl[i,:] is the Verlet list of atom i.
        * vl[i,0] = ni is the number of atoms in Verlet list of atom i.
        * vl[i,1:ni+1] contains the indices of the ni neighbouring atoms.

        This 2D data structure is linearised to facilitate passing the same
        linear data structure to Fortran and C++ in an efficient way. The
        linearised data structure consist of:

        *   vl_list : 1D numpy array containing all Verlet lists, one after the other,
                      i.e. vl(0), vl(1), ..., vl(natoms-1), with natoms the total
                      number of atoms.
        *   cl_size : 1D numpy array containing the number of atoms in the individual
                      Verlet lists, i.e. length of each Verlet list.
        *   cl_offset : 1D numpy array containing the starting position of all the
                      Verlet lists in the cl_list array.
        """
        self.cutoff = cutoff
        self.vl2d = None
        self.debug = False

    @property
    def natoms(self):
        if self.linearised():
            return len(self.vl_size)
        else:
            return len(self.vl2d)

    def allocate_2d(self, natoms):
        """allocate and initialize the 2D data structure.

        :param natoms: number of atoms to accomodate

        If the datastructure was already allocate, and large enough
        to accommodate the atoms and the contacts, it is only reinitialized.
        """

        # the maximum number of elements in a Verlet list cannot exceed natoms -1
        nneighbours = min(VerletList.nneighbours, natoms - 1)

        if hasattr(self, 'vl2d') and not self.vl2d is None:
            # reuse it
            len_vl2d = len(self.vl2d)
            if len_vl2d < natoms:
                # reset the neighbour counters
                for i in range(len_vl2d):
                    self.vl2d[i][0] = 0
                # too small, add elements at the end
                self.vl2d.extend((natoms - len_vl2d) * [[0, None]])
            elif len_vl2d > natoms:
                # too large, shrink
                self.vl2d = self.vl2d[:natoms]
                # reset the neighbour counters
                for i in range(natoms):
                    self.vl2d[i][0] = 0
        else:
            # create from scratch
            #   self.vl2d[i][0] number of neighbours of atom i
            #   self.vl2d[i][1] numpy array with the neighbour list. Initially None, but allocated
            #       when needed and grows dynamically.

            self.vl2d = natoms * [[]]
            for i in range(natoms):
                self.vl2d[i] = [0, None]

        # garbage collection of the linear data structure
        self.vl_list = None
        self.vl_size = None
        self.vl_offset = None

    def __str__(self):
        s = "verlet lists:\n"
        max_nneighbours = 0
        for i in range(self.natoms):
            vli = self.verlet_list(i)
            s += f'({i}): {vli}\n'
            max_nneighours = max(len(vli), max_nneighbours)

        s += f'max elements: {max_nneighbours}, linearised={self.linearised()}\n'
        return s

    def add(self, i, j):
        """Add pair (i,j) to the Verlet list.

        (for the 2D vl).
        """

        # increase the count of atom i
        n = self.vl2d[i][0]
        self.vl2d[i][0] += 1

        # Make sure that the Verlet list of atom i can accommodate the extra neighbour:
        if self.vl2d[i][1] is None:
            # allocate initial array
            self.vl2d[i][1] = np.empty(VerletList.nneighbours, dtype=int)
        if n == self.vl2d[i][1].shape[0]:
            # grow current array
            self.vl2d[i][1] = np.append(self.vl2d[i][1], np.empty(VerletList.nneighbours, dtype=int))

        # add the neighbour:
        self.vl2d[i][1][n] = j

    def build(self, r, keep2d=False):
        """Build the Verlet list from the positions.

        Brute force approach, but using array arithmetic, rather than pairwise
        computations. This algorithm has complexity O(N), but is significantly
        faster than ``build_simple()``.

        :param list r: list of numpy arrays with atom coordinates: r = [x, y, z]
        :param bool keep2d: if True the 2D Verlet list data structure is not deleted after linearisation.
        """
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        self.allocate_2d(len(x))
        rc2 = self.cutoff ** 2

        ri2 = np.empty((self.natoms,), dtype=r.dtype)
        rij = np.empty_like(r)
        for i in range(self.natoms - 1):
            rij[i + 1:, :] = r[i + 1:, :] - r[i, :]
            if self.debug:
                ri2 = 0
            ri2[i + 1:] = np.einsum('ij,ij->i', rij[i + 1:, :], rij[i + 1:, :])
            for j in range(i + 1, self.natoms):
                if ri2[j] <= rc2:
                    self.add(i, j)
        self.linearise(keep2d)

    def build_simple(self, r, keep2d=False):
        """Build the Verlet list from the positions.

        Brute force approach, in the simplest way.
        This algorithm has complexity O(N).

        :param list r: numpy array with atom coordinates: r.shape = (n,3)
        """
        self.allocate_2d(r.shape[0])
        rc2 = self.cutoff ** 2
        for i in range(self.natoms - 1):
            ri = r[i, :]
            for j in range(i + 1, self.natoms):
                rj = r[j, :]
                rij = rj - ri
                rij2 = np.dot(rij, rij)
                if rij2 <= rc2:
                    self.add(i, j)
        self.linearise(keep2d=keep2d)

    def build_grid(self, r, grid, keep2d=False):
        """Build Verlet lists using a grid.

        This algorithm has complexity O(N).

        :param list r: list of numpy arrays with atom coordinates: r = [x, y, z]
       """
        x = r[0]
        y = r[1]
        z = r[2]
        if not grid.linearised():
            raise ValueError("The grid list must be built and linearised first.")

        self.allocate_2d(len(x))
        rc2 = self.cutoff ** 2
        # loop over all cells
        for m in range(grid.K[2]):
            for l in range(grid.K[1]):
                for k in range(grid.K[0]):
                    cklm = grid.cell_list(k, l, m)
                    natoms_in_cklm = len(cklm)
                    # loop over all atom pairs in cklm
                    for ia in range(natoms_in_cklm):
                        i = cklm[ia]
                        for j in cklm[ia + 1:]:
                            rij2 = (x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2
                            if rij2 <= rc2:
                                self.add(i, j)
                    # loop over neighbouring cells. If the cell does not exist an IndexError is raised
                    for klm2 in ((k + 1, l, m)  # one ahead in the x-direction
                                 , (k - 1, l + 1, m)  # three ahead in the y-direction
                                 , (k, l + 1, m)
                                 , (k + 1, l + 1, m)
                                 , (k - 1, l - 1, m + 1)  # nine ahead in the z-direction
                                 , (k, l - 1, m + 1)
                                 , (k + 1, l - 1, m + 1)
                                 , (k - 1, l, m + 1)
                                 , (k, l, m + 1)
                                 , (k + 1, l, m + 1)
                                 , (k - 1, l + 1, m + 1)
                                 , (k, l + 1, m + 1)
                                 , (k + 1, l + 1, m + 1)
                                 ):
                        try:
                            cklm2 = grid.cell_list(*klm2)
                        except IndexError:
                            pass  # Cell kl2 does not exist
                        else:  # The else clause is executed only when the try clause does not raise an error
                            # loop over all atom pairs i,j with i in cklm and j in cklm2
                            for i in cklm:
                                for j in cklm2:
                                    rij2 = (x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2
                                    if rij2 <= rc2:
                                        self.add(i, j)

        self.linearise(keep2d=keep2d)

    def linearised(self):
        return not self.vl_list is None

    def linearise(self, keep2d=False):
        """linearise the Verlet list.

        :param bool keep2d: keep self.vl2d or not. True is used for testing purposes.
        """
        natoms = self.natoms
        nneighbours_total = 0
        for i in range(natoms):
            nneighbours_total += self.vl2d[i][0]

        self.vl_size = np.empty(natoms, dtype=int)
        self.vl_offset = np.empty(natoms, dtype=int)
        self.vl_list = np.empty(nneighbours_total, dtype=int)

        offset = 0
        for i in range(natoms):
            nneighbours_i = self.vl2d[i][0]
            self.vl_size[i] = nneighbours_i
            self.vl_offset[i] = offset
            if nneighbours_i:
                self.vl_list[offset:offset + nneighbours_i] = self.vl2d[i][1][:nneighbours_i]
            offset += nneighbours_i

        if not keep2d:
            self.vl2d = None  # garbage collection takes care of it.

    def verlet_list(self, i):
        """Return the Verlet list of atom i.

        :return: view of a numpy array
        """
        if self.linearised():
            offset = self.vl_offset[i]
            nneighbours = self.vl_size[i]
            vli = self.vl_list[offset:offset + nneighbours]
        else:
            nneighbours = self.vl2d[i][0]
            vli = self.vl2d[i,][0][:nneighbours]

        return vli

    def has(self, ij):
        """Test if the verlet list of atom i contains atom j, ij = (i,j).

        :param tuple ij: (i,j)
        """
        i, j = ij
        return j in self.verlet_list(i)

def vl2set(vl):
    """Convert VerletList object into set of pairs."""
    pairs = set()
    if isinstance(vl, VerletList):
        for i in range(vl.natoms):
            vli = vl.verlet_list(i)
            n_pairs_i = len(vli)
            for j in vli:
                pair = (i, j) if i < j else (j, i)
                # print(pair)
                pairs.add(pair)
    else:
        # C++ implementation of Verlet list: VList
        for i in range(vl.natoms()):
            for k in range(vl.ncontacts(i)):
                j = vl.contact(i, k)
                pair = (i, j) if i < j else (j, i)
                # print(pair)
                pairs.add(pair)

    return pairs
