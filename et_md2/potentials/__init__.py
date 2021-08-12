# -*- coding: utf-8 -*-

"""
Module et_md2.potentials
========================

A submodule for interaction potentials.

"""

class LJ_py:
    """Lennard-Jones potentiol

    See https://en.wikipedia.org/wiki/Lennard-Jones_potential.

    V_lj(r) = 4*epsilon*[ (sigma/r)**12 - (sigma/r)**6 ]

    :param float epsilon: epsilon
    :param float sigma: sigma
    """
    def __init__(self,epsilon=0.25, sigma=1):
        self.four_epsilon = 4.0*epsilon
        self.inv_sigma = (1.0/sigma)


    def r0(self):
        """Equilibrium distance."""
        return pow(2.,1/6) / self.inv_sigma


    def interaction_energy(self, rij2):
        """Compute the Lennard-Jones potential

        :param float|np.array rij2: squared distance between atoms.
        :returns: a float.
        """
        rij_sigma_2 = rij2 * self.inv_sigma * self.inv_sigma
        rm6 = 1./(rij_sigma_2*rij_sigma_2*rij_sigma_2)
        vlj = self.four_epsilon*(rm6 - 1.0)*rm6
        return vlj


    def force_factor(self, rij2):
        """Lennard-Jones force magnitude exerted by atom j on atom i.

        :param float|np.array rij2: squared interatomic distance from atom i to atom j
        :return: fij
        """
        rij_sigma_2 = rij2 * self.inv_sigma * self.inv_sigma
        rm2 = 1.0 / rij_sigma_2
        rm6 = (rm2 * rm2 * rm2)
        f = self.four_epsilon * self.inv_sigma * 6.0 * (1.0 - 2.0 * rm6) * rm6 * rm2
        return f


    def force(self, rij):
        """Lennard-Jones force exerted by atom j on atom i.

        :param 3-tuple of float|np.array rij: vector from atom i to atom j
        :return: Fij, list of float|np.array
        """
        rij2 = rij[0] ** 2 + rij[1] ** 2 + rij[2] ** 2
        f = self.force_factor(rij2)
        return ( f * rij[0]
               , f * rij[1]
               , f * rij[2]
               )
