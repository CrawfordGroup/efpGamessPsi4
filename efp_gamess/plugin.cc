/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. 
 */

#include <memory>

#include <psi4/psi4-dec.h>

#include <psi4/libmints/basisset.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libpsio/psio.hpp>

namespace psi {
namespace efp_gamess {

extern "C" PSI_API int read_options(std::string name, Options &options) {
  if (name == "EFP_GAMESS" || options.read_globals()) {
    /*- The amount of information printed to the output file -*/
    options.add_int("PRINT", 1);
    options.add_str("TRANS_MAT", "");
    options.add_str_i("HDF5_FILENAME", "");
  }

  return true;
}

extern "C" PSI_API SharedMatrix efp_gamess(SharedWavefunction ref_wfn,
                                           Options &options) {
  int print = options.get_int("PRINT");
  std::string which_trans = options.get_str("TRANS_MAT");
  auto basis = ref_wfn->basisset();
  int nao = basis->nao();
  int nso = ref_wfn->nso();
  int nmo = ref_wfn->nmo();

  auto UC_ = std::make_shared<Matrix>("C Transformation matrix", nso, nao);
  auto UH_ = std::make_shared<Matrix>("H Transformation matrix", nso, nao);
  bool puream = basis->has_puream();
  int so_off = 0;
  int ao_off = 0;
  for (int shell = 0; shell < basis->nshell(); ++shell) {
    int am = basis->shell(shell).am();
    /*
     * Re-map the SO indices.  The Cartesian ordering is as follows.
     *
     *         | 0 |   1   |         2         |                   3
     * --------|---|-------|-------------------|----------------------------------------
     *  Psi4:  | s | x y z | xx xy xz yy yz zz | xxx xxy xxz xyy xyz xzz yyy yyz
     * yzz zzz GAMESS | s | x y z | xx yy zz xy xz yz | xxx yyy zzz xxy xxz xyy
     * yyz xzz yzz xyz
     * --------|---|-------|-------------------|----------------------------------------
     *  offset | 0 | 0 1 2 |  0  3  4  1  5  2 |  0   3   4   5   9   7   1   6
     * 8   2
     *
     * We also need to account for the missing normalization factor
     *      __________________________
     *     /         2^(2l)
     *    / --------------------------
     *  \/ (2lx-1)!!(2ly-1)!!(2lz-1)!!
     *
     * For sphericals, we need to do a little more.  It looks like GAMESS always
     * reports spherical harmonics in terms of their Cartesian components, so we
     * just need to do a quick remapping, as follows.  The Psi4 ordering is
     * always 1, 1c, 1s, 2c, 2s,...
     *
     * l=1  | x | y | z
     * -----------------
     * 1_0  | 0 | 0 | 1
     * 1_1c | 1 | 0 | 0
     * 1_1s | 0 | 1 | 0
     *
     * l=2  | xx        | xy     | yy         | xz      | yz     | zz
     * --------------------------------------------------------------
     * 2_0  | -1/2      | 0      | -1/2       | 0       | 0      | 1
     * 2_1c | 0         | 0      | 0          | sqrt(3) | 0      | 0
     * 2_1s | 0         | 0      | 0          | 0       | sqrt(3)| 0
     * 2_2c | sqrt(3)/2 | 0      | -sqrt(3)/2 | 0       | 0      | 0
     * 2_2s | 0         | sqrt(3)| 0          | 0       | 0      | 0
     *
     * l=3 | xxx | xxy | xyy | yyy | xxz | xyz | yyz | xzz | yzz | zzz
     * ----------------------------------------------------------------
     * 3_0 | 0 | 0 | 0 | 0 | -frac32 | 0 | -frac32 | 0 | 0 | 1
     * 3_1c| -sqrtfrac322 | 0 | -fracsqrtfrac322 | 0 | 0 | 0 | 0 | sqrt6 | 0 | 0
     * 3_1s| 0 | -fracsqrtfrac322 | 0 | -fracsqrtfrac322 | 0 | 0 | 0 | 0 | sqrt6
     * | 0 3_2c| 0 | 0 | 0 | 0 | fracsqrt152 | 0 | -fracsqrt152 | 0 | 0 | 0
     * 3_2s| 0 | 0 | 0 | 0 | 0 | sqrt15 | 0 | 0 | 0 | 0
     * 3_3c| fracsqrtfrac522 | 0 | -frac3 sqrtfrac522 | 0 | 0 | 0 | 0 | 0 | 0 |
     * 0 3_3s| 0 | frac3 sqrtfrac522 | 0 | -fracsqrtfrac522 | 0 | 0 | 0 | 0 | 0
     * | 0
     */
    if (puream) {
      if (am == 0) {
        // s
        UC_->set(so_off, ao_off, 1.0);
        UH_->set(so_off, ao_off, 1.0);
        so_off++;
        ao_off++;
      } else if (am == 1) {
        // 1_0  <- z
        UC_->set(so_off + 0, ao_off + 2, 1.0);
        UH_->set(so_off + 0, ao_off + 2, 1.0);
        // 1_1c  <- x
        UC_->set(so_off + 1, ao_off + 0, 1.0);
        UH_->set(so_off + 1, ao_off + 0, 1.0);
        // 1_1s  <- y
        UC_->set(so_off + 2, ao_off + 1, 1.0);
        UH_->set(so_off + 2, ao_off + 1, 1.0);
        so_off += 3;
        ao_off += 3;
      } else if (am == 2) {
        double Caa = 1.0;
        double Cab = sqrt(3.0);
        double Haa = 1.0;
        double Hab = 1.0 / sqrt(3.0);

        double C0 = 2.0 / 3.0;
        double C1c = 1.0 / sqrt(3.0);
        double C1s = 1.0 / sqrt(3.0);
        double C2c = 1.0 / sqrt(3.0);
        double C2s = 1.0 / sqrt(3.0);
        double H0 = 1.0;
        double H1c = sqrt(3.0);
        double H1s = sqrt(3.0);
        double H2c = sqrt(3.0) / 2.0;
        double H2s = sqrt(3.0);
        // 2_0 <- zz - 0.5xx - 0.5yy
        UC_->set(so_off + 0, ao_off + 2, 1.0 * C0 * Caa);
        UC_->set(so_off + 0, ao_off + 0, -0.5 * C0 * Caa);
        UC_->set(so_off + 0, ao_off + 1, -0.5 * C0 * Caa);
        UH_->set(so_off + 0, ao_off + 2, 1.0 * H0 * Haa);
        UH_->set(so_off + 0, ao_off + 0, -0.5 * H0 * Haa);
        UH_->set(so_off + 0, ao_off + 1, -0.5 * H0 * Haa);
        // 2_1c <- xz
        UC_->set(so_off + 1, ao_off + 4, 1.0 * C1c * Cab);
        UH_->set(so_off + 1, ao_off + 4, 1.0 * H1c * Hab);
        // 2_1s <- yz
        UC_->set(so_off + 2, ao_off + 5, 1.0 * C1s * Cab);
        UH_->set(so_off + 2, ao_off + 5, 1.0 * H1s * Hab);
        // 2_2c <- xx - yy
        UC_->set(so_off + 3, ao_off + 0, 1.0 * C2c * Haa);
        UC_->set(so_off + 3, ao_off + 1, -1.0 * C2c * Haa);
        UH_->set(so_off + 3, ao_off + 0, 1.0 * H2c * Haa);
        UH_->set(so_off + 3, ao_off + 1, -1.0 * H2c * Haa);
        // 2_2s <- xy
        UC_->set(so_off + 4, ao_off + 3, 1.0 * C2s * Cab);
        UH_->set(so_off + 4, ao_off + 3, 1.0 * H2s * Hab);
        so_off += 5;
        ao_off += 6;
      } else if (am == 3) {
        double Caaa = 1.0;
        double Caab = sqrt(5.0);
        double Cabc = sqrt(15.0);
        double Haaa = 1.0;
        double Haab = 1.0 / sqrt(5.0);
        double Habc = 1.0 / sqrt(15.0);

        double C0 = 2.0 / 11.0;        // good
        double C1c = 1.0 / sqrt(15.0); // bad
        double C1s = 1.0 / sqrt(15.0); // bad
        double C2c = 1.0 / sqrt(15.0); // good
        double C2s = 1.0 / sqrt(15.0); // good
        double C3c = 1.0 / sqrt(15.0); // bad
        double C3s = 1.0 / sqrt(15.0); // bad
        double H0 = 1.0;
        double H1c = sqrt(6.0);
        double H1s = sqrt(6.0);
        double H2c = sqrt(15.0 / 4.0);
        double H2s = sqrt(15.0);
        double H3c = sqrt(5.0 / 2.0);
        double H3s = sqrt(5.0 / 2.0);
        //   xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
        //    0   1   2   3   4   5   6   7   8   9
        // 3_0 <- zzz - 1.5xxz - 1.5yyz
        UC_->set(so_off + 0, ao_off + 2, 1.0 * C0 * Caaa);
        UC_->set(so_off + 0, ao_off + 4, -1.5 * C0 * Caab);
        UC_->set(so_off + 0, ao_off + 6, -1.5 * C0 * Caab);
        UH_->set(so_off + 0, ao_off + 2, 1.0 * H0 * Haaa);
        UH_->set(so_off + 0, ao_off + 4, -1.5 * H0 * Haab);
        UH_->set(so_off + 0, ao_off + 6, -1.5 * H0 * Haab);
        // 3_1c <- sqrt(6)xzz - sqrt(3/8)xyy - sqrt(3/8)xxx
        UC_->set(so_off + 1, ao_off + 7, 1.0 * C1c * Caab);
        UC_->set(so_off + 1, ao_off + 5, -0.25 * C1c * Caab);
        UC_->set(so_off + 1, ao_off + 0, -0.25 * C1c * Caaa);
        UH_->set(so_off + 1, ao_off + 7, 1.0 * H1c * Haab);
        UH_->set(so_off + 1, ao_off + 5, -0.25 * H1c * Haab);
        UH_->set(so_off + 1, ao_off + 0, -0.25 * H1c * Haaa);
        // 3_1s <- sqrt(6)yzz - sqrt(3/8)yyy - sqrt(3/8)xxy
        UC_->set(so_off + 2, ao_off + 8, 1.0 * C1s * Caab);
        UC_->set(so_off + 2, ao_off + 1, -0.25 * C1s * Caaa);
        UC_->set(so_off + 2, ao_off + 3, -0.25 * C1s * Caab);
        UH_->set(so_off + 2, ao_off + 8, 1.0 * H1s * Haab);
        UH_->set(so_off + 2, ao_off + 1, -0.25 * H1s * Haaa);
        UH_->set(so_off + 2, ao_off + 3, -0.25 * H1s * Haab);
        // 3_2c <- sqrt(15/4)xxz - sqrt(15/4)yyz
        UC_->set(so_off + 3, ao_off + 4, 1.0 * C2c * Caab);
        UC_->set(so_off + 3, ao_off + 6, -1.0 * C2c * Caab);
        UH_->set(so_off + 3, ao_off + 4, 1.0 * H2c * Haab);
        UH_->set(so_off + 3, ao_off + 6, -1.0 * H2c * Haab);
        // 3_2s <- sqrt(15)xyz
        UC_->set(so_off + 4, ao_off + 9, 1.0 * C2s * Cabc);
        UH_->set(so_off + 4, ao_off + 9, 1.0 * H2s * Habc);
        // 3_3c <- sqrt(5/2)/2xxx - 3sqrt(5/2)/2xyy
        UC_->set(so_off + 5, ao_off + 0, 0.5 * C3c * Caaa);
        UC_->set(so_off + 5, ao_off + 5, -1.5 * C3c * Caab);
        UH_->set(so_off + 5, ao_off + 0, 0.5 * H3c * Haaa);
        UH_->set(so_off + 5, ao_off + 5, -1.5 * H3c * Haab);
        // 3_3s <- 3sqrt(5/2)/2xxy - sqrt(5/2)/2yyy
        UC_->set(so_off + 6, ao_off + 3, 1.5 * C3s * Caab);
        UC_->set(so_off + 6, ao_off + 1, -0.5 * C3s * Caaa);
        UH_->set(so_off + 6, ao_off + 3, 1.5 * H3s * Haab);
        UH_->set(so_off + 6, ao_off + 1, -0.5 * H3s * Haaa);
        so_off += 7;
        ao_off += 10;
      } else {
        throw PSIEXCEPTION("f functions not yet implemented for pure A.M.");
      }
    } else {
      if (am == 0) {
        // s
        UC_->set(so_off, ao_off, 1.0);
        UH_->set(so_off, ao_off, 1.0);
        so_off++;
        ao_off++;
      } else if (am == 1) {
        // x
        UC_->set(so_off + 0, ao_off + 0, 1.0);
        UH_->set(so_off + 0, ao_off + 0, 1.0);
        // y
        UC_->set(so_off + 1, ao_off + 1, 1.0);
        UH_->set(so_off + 1, ao_off + 1, 1.0);
        // z
        UC_->set(so_off + 2, ao_off + 2, 1.0);
        UH_->set(so_off + 2, ao_off + 2, 1.0);
        so_off += 3;
        ao_off += 3;
      } else if (am == 2) {
        double Caa = 1.0;
        double Cab = sqrt(3.0);
        double Haa = 1.0;
        double Hab = 1.0 / sqrt(3.0);
        // xx
        UC_->set(so_off + 0, ao_off + 0, Caa);
        UH_->set(so_off + 0, ao_off + 0, Haa);
        // xy
        UC_->set(so_off + 1, ao_off + 3, Cab);
        UH_->set(so_off + 1, ao_off + 3, Hab);
        // xz
        UC_->set(so_off + 2, ao_off + 4, Cab);
        UH_->set(so_off + 2, ao_off + 4, Hab);
        // yy
        UC_->set(so_off + 3, ao_off + 1, Caa);
        UH_->set(so_off + 3, ao_off + 1, Haa);
        // yz
        UC_->set(so_off + 4, ao_off + 5, Cab);
        UH_->set(so_off + 4, ao_off + 5, Hab);
        // zz
        UC_->set(so_off + 5, ao_off + 2, Caa);
        UH_->set(so_off + 5, ao_off + 2, Haa);
        so_off += 6;
        ao_off += 6;
      } else if (am == 3) {
        double Caaa = 1.0;
        double Caab = sqrt(5.0);
        double Cabc = sqrt(15.0);
        double Haaa = 1.0;
        double Haab = 1.0 / sqrt(5.0);
        double Habc = 1.0 / sqrt(15.0);
        // xxx
        UC_->set(so_off + 0, ao_off + 0, Caaa);
        UH_->set(so_off + 0, ao_off + 0, Haaa);
        // xxy
        UC_->set(so_off + 1, ao_off + 3, Caab);
        UH_->set(so_off + 1, ao_off + 3, Haab);
        // xxz
        UC_->set(so_off + 2, ao_off + 4, Caab);
        UH_->set(so_off + 2, ao_off + 4, Haab);
        // xyy
        UC_->set(so_off + 3, ao_off + 5, Caab);
        UH_->set(so_off + 3, ao_off + 5, Haab);
        // xyz
        UC_->set(so_off + 4, ao_off + 9, Cabc);
        UH_->set(so_off + 4, ao_off + 9, Habc);
        // xzz
        UC_->set(so_off + 5, ao_off + 7, Caab);
        UH_->set(so_off + 5, ao_off + 7, Haab);
        // yyy
        UC_->set(so_off + 6, ao_off + 1, Caaa);
        UH_->set(so_off + 6, ao_off + 1, Haaa);
        // yyz
        UC_->set(so_off + 7, ao_off + 6, Caab);
        UH_->set(so_off + 7, ao_off + 6, Haab);
        // yzz
        UC_->set(so_off + 8, ao_off + 8, Caab);
        UH_->set(so_off + 8, ao_off + 8, Haab);
        // zzz
        UC_->set(so_off + 9, ao_off + 2, Caaa);
        UH_->set(so_off + 9, ao_off + 2, Haaa);
        so_off += 10;
        ao_off += 10;
      } else {
        throw PSIEXCEPTION("f functions not yet implemented for pure A.M.");
      }
    }
  }

  // SharedMatrix diffMat = SharedMatrix(new Matrix("Difference in trans.
  // mats",nso,nao)); diffMat->add(UC_); diffMat->subtract(UH_);
  // diffMat->print_out();

  if (which_trans == "C") {
    outfile->Printf("Coleman I am return UC_\n");
    return UC_;
  } else if (which_trans == "F") {
    outfile->Printf("Coleman I am return UH_\n");
    return UH_;
  } else {
    outfile->Printf("Coleman, I got an invalid value for option trans_mat\n");
    exit(1);
  }
}

} // namespace efp_gamess
} // namespace psi
