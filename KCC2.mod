:------------------------------------------------------------------------------
TITLE Model of KCC2 regulation
:------------------------------------------------------------------------------

NEURON {
  SUFFIX KCC2
  USEION ca READ cai
  USEION mkcc2 WRITE mkcc2i VALENCE 1
  RANGE A_M, B_M
  RANGE R_M, R_MP
  RANGE R_K, R_P, V_K, V_P, H_K, H_P, B_K, B_P
  RANGE kin_active, kin_inactive, cyt, memb, membp, v_K, v_P
}

UNITS {
  (M)   = (1/liter)
  (mM)  = (milliM)
}

PARAMETER {
  A_M  = 1  (/s)    < 0, 1e9 >    : rate of transfer from cytosolic to membrane-bound KCC2
  B_M  = 1  (/s)    < 0, 1e9 >    : rate of transfer from membrane-bound to cytosolic KCC2
  R_M  = 5  (/s)    < 0, 1e9 >    : KCC2 phosphorylation coefficient
  R_MP = 5  (/s)    < 0, 1e9 >    : KCC2 dephosphorylation coefficient
  R_K  = .5 (mM)    < 0, 1e9 >    : concent. of Ca2+ at which kinase activation rate is half max
  R_P  = .5 (mM)    < 0, 1e9 >    : concent. of Ca2+ at which phosphotase activation rate is half max
  V_K  = 5  (/s)    < 0, 1e9 >    : maximum rate of kinase activation
  V_P  = 5  (/s)    < 0, 1e9 >    : maximum rate of phosphotase activation
  H_K  = 2  (1)     < 0, 1e9 >    : Hill coefficient for kinase activation equation
  H_P  = 2  (1)     < 0, 1e9 >    : Hill coefficient for phosphotase activation equation
  B_K  = 5  (/s)    < 0, 1e9 >    : rate of kinase inactivation
  B_P  = 5  (/s)    < 0, 1e9 >    : rate of phosphotase inactivation
}

ASSIGNED {
  mkcc2i (mM) : concentration of membrane-bound KCC2
  cai    (mM) : calcium concentration
  v_K    (/s) : rate of kinase activation
  v_P    (/s) : rate of phosphotase activation
}

STATE {
  kin_active
  kin_inactive
  phos_active
  phos_inactive
  cyt
  memb
  membp
}

INITIAL {
  kin_active    = .00004
  kin_inactive  = .99996
  phos_active   = .00004
  phos_inactive = .99996
  cyt           = .891
  memb          = .05
  membp         = .059
}

BREAKPOINT {
  SOLVE states METHOD cnexp
}

DERIVATIVE states {
  enzymes(cai)

  cyt'   = (1e-3)*(B_M*memb - A_M*cyt)
  memb'  = (1e-3)*(A_M*cyt - (B_M + R_MP*kin_active)*memb + R_M*phos_active*membp)
  membp' = (1e-3)*(-R_M*phos_active*membp + R_MP*kin_active*memb)

  kin_inactive   = 1 - kin_active
  phos_inactive  = 1 - phos_active
  kin_active'    = (1e-3)*(-R_MP*memb*kin_active + kin_inactive*v_K - B_K*kin_active)
  phos_active'   = (1e-3)*(-R_M*membp*phos_active + phos_inactive*v_P - B_P*phos_active)

  if (memb + membp > 0) {
    mkcc2i = (memb + membp)*1(mM)
  } else {
    mkcc2i = 0
  }
}

PROCEDURE enzymes (Caint (mM)) {
  v_K = (((Caint/1(mM))^H_K)*V_K)/((R_K/1(mM))^H_K + (Caint/1(mM))^H_K)
  v_P = (((Caint/1(mM))^H_P)*V_P)/((R_P/1(mM))^H_P + (Caint/1(mM))^H_P)
}