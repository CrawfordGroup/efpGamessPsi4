#
# @BEGIN LICENSE
#
# gamess_h5 by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np

import h5py
import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util


def run_efp_gamess(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    efp_gamess can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('efp_gamess')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('EFP_GAMESS', 'PRINT', 1)
    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    # This prints a bound python method

    # Get C transformation matrix from plugin
    psi4.core.set_local_option('EFP_GAMESS', 'TRANS_MAT', 'C')
    trans_mat_c = psi4.core.plugin('efp_gamess.so', ref_wfn)
    #print("\nPrinting C transformation matrix (as array)")
    #print(np.asarray(trans_mat_c))
    #trans_mat_c.print_out()

    # Get F transformation matrix from plugin
    psi4.core.set_local_option('EFP_GAMESS', 'TRANS_MAT', 'F')
    trans_mat_f = psi4.core.plugin('efp_gamess.so', ref_wfn)
    #print("\nPrinting F transformation matrix (as array)")
    #print(np.asarray(trans_mat_f))
    #trans_mat_f.print_out()

    # Define hdf5 file
    f = h5py.File("form.h5", "r")
    # Group is "EFPcalc"
    group = f["EFPcalc"]
    #print("\nListing dataset in h5 file EFPcalc group")
    #print(list(group.keys()))
    # Define dataset for Fock matrix
    fock_dset = group['CONVERGED TOTAL FOCK MATRIX']
    fock_np_lt = np.array(fock_dset)
    #print("\nGAMESS Fock matrix as numpy array (lower triangle)")
    #print(fock_np_lt)
    nbf = ref_wfn.basisset().nbf()
    nao = ref_wfn.basisset().nao()
    nmo = ref_wfn.nmo()
    nso = ref_wfn.nso()
    # Print basis set info for debugging for now
    #print("Basis set info")
    #print(ref_wfn.basisset().has_puream())
    # Make iterator to turn lower triangle Fock from Gamess into full matrix
    fock_lt_iter = np.nditer(fock_np_lt, order='C')
    #print("nbf = ",nbf)
    fock_np = np.zeros((nao, nao))
    fock_np.shape = (nao, nao)
    for i in range(nao):
        for j in range(i+1):
            fock_np[i][j] = fock_lt_iter[0]
            fock_np[j][i] = fock_lt_iter[0]
            fock_lt_iter.iternext()

    #print("Upper triangle Fock from Gamess turned into full Matrix:")
    #print(fock_np)
    # Define data set for MO coefficients
    # Mo dset array is (146,154), F is (154,154)
    # num ao is 154, # bf is 146
    # Cgamess = (nao,nmo)
    mo_dset = group['MO_coeff']
    mo_np = np.transpose(np.array(mo_dset))
    #mo_np=np.array(mo_dset)
    #print("\nGAMESS MO coeficients")
    #print(mo_np)

    # Transform Gamess MO ordering to PSI4 order
    #print("\nPSI4 MO coefficients")
    psi4.core.print_out("Transforming MO coefficients\n")
    psi4_C = np.einsum("ik,kj->ij", trans_mat_c, mo_np)
    psi4.core.print_out("DONE transforming MO coefficients\n")
    #print("\nTransformed MO coefficients")
    #print(psi4_C)

    #TODO transpose F when getting from Gamess fortran
    # Transform Gamess Fock matrix ordering to PSI4 order
    psi4.core.print_out("Transforming GAMESS Fock matrix\n")
    psi4_F = np.einsum("ik,kl,jl->ij", trans_mat_f, fock_np, trans_mat_f)
    psi4.core.print_out("DONE transforming GAMESS Fock matrix\n")
    #psi4_F_tmp = np.matmul(fock_np,np.transpose(trans_mat_f))
    #psi4_F = np.matmul(trans_mat_f, psi4_F_tmp)
    #print(psi4_F)

    # This returns psi4.core.Matrix objects
    Fa = ref_wfn.Fa()
    Fa.copy(psi4.core.Matrix.from_array(psi4_F))
    Ca = ref_wfn.Ca()
    Ca.copy(psi4.core.Matrix.from_array(psi4_C))
    #ref_wfn.Fa().print_out()
    #ref_wfn.Ca().print_out()

    return ref_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['efp_gamess'] = run_efp_gamess
