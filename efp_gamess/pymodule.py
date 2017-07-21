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

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
import h5py
import numpy as np

def run_efp_gamess(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    efp_gamess can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('efp_gamess')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)

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
    print("\nPrinting ref_wfn.Fa")
    print(ref_wfn.Fa)
    # This defines a psi4.core.Matrix object and prints it
    print(np.asarray(ref_wfn.Fa()))
    #Farray=np.asarray(ref_wfn.Fa())
    Farray=ref_wfn.Fa().to_array()
    print("\nPrinting Farray, which is np.asarray(ref_wfn.Fa())")
    print(Farray)

    print("\nPrinting Farray asarray")
    print(np.asarray(Farray))
    # I think this does nothing different from the above
    Fa=ref_wfn.Fa()
    print("\nPrinting Fa=ref_wfn.Fa()")
    print(Fa)
    print("\nPrinting np.asarray(Fa)")
    print(np.asarray(Fa))
    # This prints to the output the matrix values
    Fa.print_out()
    # This actually changes Fa and the ref_wfn Fa value
    Fa.set(0,0,0.0)
    Fa.print_out()
    #ref_wfn.Fa=Fa
    # Tests my previous assertion about ref_wfn changes
    Fa_zeroed=ref_wfn.Fa()
    Fa_zeroed.print_out()

    # Change return name to the trans matrix I'm trying to send back
    #efp_gamess_wfn = psi4.core.plugin('efp_gamess.so', ref_wfn)

    
    psi4.core.set_local_option('efp_gamess', 'coleman', 0)

    trans_mat_c=psi4.core.plugin('efp_gamess.so',ref_wfn)
    print("\nPrinting C transformation matrix")
    print(trans_mat_c)
    print("\nPrinting C transformation matrix (as array)")
    print(np.asarray(trans_mat_c))
    trans_mat_c.print_out()

    psi4.core.print_out("\nalright coleman H coming now")
    psi4.core.set_local_option('efp_gamess', 'coleman', 1)

    trans_mat_h=psi4.core.plugin('efp_gamess.so',ref_wfn)
    print("\nPrinting H transformation matrix")
    print(trans_mat_h)
    print("\nPrinting H transformation matrix (as array)")
    print(np.asarray(trans_mat_h))
    trans_mat_h.print_out()

    # Test adding these together
    print("Adding np.asarray(trans_mat_h)+np.asarray(trans_mat_c)\n")
    test_sum=np.asarray(trans_mat_h)+np.asarray(trans_mat_c)
    print(test_sum)
    #test_sum2=np.asarray(Fa)+np.asarray(trans_mat_c)
    #print(test_sum2)
    #test_sum.print_out()
    # Find out what trans_mat is (it's like the others)
    #print(trans_mat)
    # print actual values
    #trans_mat.print_out()

    # Change return to original ref_wfn, hopefully changed
    #return efp_gamess_wfn
    return ref_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['efp_gamess'] = run_efp_gamess


def exampleFN():
    # Your Python code goes here
    pass
