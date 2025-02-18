/*
Copyright (c) 1999-2004 I. B. Smirnov

The file can be used, copied, modified, and distributed
according to the terms of GNU Lesser General Public License version 2.1
as published by the Free Software Foundation,
and provided that the above copyright notice, this permission notice,
and notices about any modifications of the original text
appear in all copies and in supporting documentation.
The file is provided "as is" without express or implied warranty.
*/
#include "wcpplib/safetl/AbsArr.h"

namespace Heed {

long max_qel_DynLinArr = 100000000;

int gconfirm_ind(const DynLinArr<long>& qel, const DynLinArr<long>& ind) {
  if (qel.get_qel() != ind.get_qel()) {
    mcerr << "gconfirm_ind(...): "
          << "qel.get_qel()!= ind.get_qel()\n"
          << "qel.get_qel()=" << qel.get_qel()
          << "ind.get_qel()=" << ind.get_qel() << '\n';
    spexit(mcerr);
  }
  long qd = qel.get_qel();
  // if( ind.get_qel() < qd) qd=ind.get_qel();
  long n;
  for (n = 0; n < qd; n++)
    if (ind[n] < 0 || ind[n] >= qel[n]) return 0;
  return 1;
}
int gconfirm_ind_ext(const DynLinArr<long>& qel, const DynLinArr<long>& ind) {
  if (qel.get_qel() > ind.get_qel()) {
    mcerr << "gconfirm_ind_ext(...): "
          << "qel.get_qel()> ind.get_qel()\n"
          << "qel.get_qel()=" << qel.get_qel()
          << " ind.get_qel()=" << ind.get_qel() << '\n';
    spexit(mcerr);
  }
  long qd = qel.get_qel();
  long n;
  for (n = 0; n < qd; n++)
    if (ind[n] < 0 || ind[n] >= qel[n]) return 0;
  for (n = qd; n < ind.get_qel(); n++)
    if (ind[n] != 0) return 0;
  return 1;
}

int find_next_comb(const DynLinArr<long>& qel, DynLinArr<long>& f) {
  long n;
  long qdim = qel.get_qel();
  if (qdim <= 0) return 0;
  if (qdim != f.get_qel()) return 0;  //@@
#ifdef DEBUG_DYNARR
  for (n = qdim - 1; n >= 0; n--) {
    if (f[n] < qel[n] - 1) {
      f[n]++;
      return 1;
    } else {
      f[n] = 0;
    }  // the first element

  }                                                  // it was last combination
  for (n = 0; n < qdim - 1; n++) f[n] = qel[n] - 1;  // the last element
  f[qdim - 1] = qel[qdim - 1];                       // next after last
#else
  for (n = qdim - 1; n >= 0; n--) {
    if (f.acu(n) < qel.acu(n) - 1) {
      f.acu(n)++;
      return 1;
    } else {
      f.acu(n) = 0;
    }  // the first element

  }  // it was last combination
  for (n = 0; n < qdim - 1; n++) f.acu(n) = qel.acu(n) - 1;  // the last element
  f.acu(qdim - 1) = qel.acu(qdim - 1);                       // next after last
#endif
  return 0;
}

int find_next_comb_not_less(const DynLinArr<long>& qel, DynLinArr<long>& f) {
  long n;
  long qdim = qel.get_qel();
  if (qdim <= 0) return 0;
  if (qdim != f.get_qel()) return 0;  //@@
  for (n = qdim - 1; n >= 0; n--) {
    if (f[n] < qel[n] - 1) {
      f[n]++;
      int n1;
      for (n1 = n + 1; n1 < qdim; n1++) f[n1] = f[n];
      return 1;
    }
  }                                                  // it was last combination
  for (n = 0; n < qdim - 1; n++) f[n] = qel[n] - 1;  // the last element
  f[qdim - 1] = qel[qdim - 1];                       // next after last
  return 0;
}

int find_prev_comb(const DynLinArr<long>& qel, DynLinArr<long>& f) {
  long n;
  long qdim = qel.get_qel();
  if (qdim <= 0) return 0;
  if (qdim != f.get_qel()) return 0;  //@@
  for (n = qdim - 1; n >= 0; n--) {
    if (f[n] >= 1) {
      f[n]--;
      return 1;
    } else {
      f[n] = qel[n] - 1;
    }  // the last element
  }
  for (n = 0; n < qdim - 1; n++) f[n] = 0;  // the first element
  f[qdim - 1] = -1;
  return 0;  // previous before first
}

DynLinArr<long> qel_communicat;
}
