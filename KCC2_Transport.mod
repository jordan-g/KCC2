NEURON {
  SUFFIX KCC2_Transport
  USEION k READ ki, ko VALENCE 1
  USEION cl READ icl, ecl, clo WRITE cli VALENCE -1
  USEION mkcc2 READ mkcc2i VALENCE 1
  RANGE axD, gtonic, itonic
  RANGE R_T, r_T
  RANGE D, v_T, V_T, transport
  GLOBAL celsius
}

UNITS {
  (mV) = (millivolt)
  (um) = (micron)
  (mA) = (milliamp)
  (M)  = (1/liter)
  (mM) = (milliM)
  (J)  = (joule)
  (S)  = (siemens)

  F   = (faraday) (coulombs)
  R   = (k-mole) (joule/degC)
  PI  = (pi) (1)
}

PARAMETER {
  axD    = 1    (um2/ms) < 0, 1e9 > : axial chloride diffusion constant
  gtonic = 0.01 (S/cm2)  < 0, 1e9 > : tonic conductance
  R_T    = 5    (mM)     < 0, 1e9 > : Michaelis constant for KCC2 ion transport
  r_T    = 5    (/s)     < 0, 1e9 > : rate of transport of KCC2-bound ions across the membrane
}

ASSIGNED {
  v         (mV)     : membrane voltage
  diam      (um)     : compartment diameter
  L         (um)     : compartment length
  icl       (mA/cm2) : chloride current
  ik        (mA/cm2) : potassium current
  itonic    (mA/cm2) : tonic GABA current
  clo       (mM)     : external chloride concentration
  ko        (mM)     : external potassium concentration
  ecl       (mV)     : chloride reversal potential
  celsius   (degC)   : temperature
  mkcc2i    (mM)     : concentration of membrane-bound KCC2
  D         (mV)     : KCC2 driving force direction
  V_T       (mM/s)   : max rate of ion transport by KCC2
  v_T       (mM/s)   : rate of ion transport by KCC2
  transport (mM/s)   : rate of ion transport accounting for direction
}

STATE {
  cli (mM) : internal chloride concentration
  ki  (mM) : internal potassium concentration
}

INITIAL {
  itonic    = 0
  D         = 0
  V_T       = 0
  v_T       = 0
  transport = 0
}

BREAKPOINT {
  SOLVE state METHOD sparse
}

KINETIC state {
  tonicgaba(v,ecl)
  transport_mm(cli, ki, mkcc2i)
  driveforce(cli,clo, ki, ko)
  
  COMPARTMENT                PI*diam*diam/4 { cli }
  LONGITUDINAL_DIFFUSION axD*PI*diam*diam/4 { cli }
  
  : synaptic chloride current - KCC2 efflux + tonic GABA current
  ~ cli << ((1e4)*icl*diam*PI/F + (1e-3)*transport*PI*diam*diam/4 + (1e4)*itonic*diam*PI/F)
}

PROCEDURE tonicgaba(Vm (mV), Em (mV)) {
  itonic = gtonic * (Vm - Em)
}

PROCEDURE transport_mm(Cint (mM), Kint (mM), Mkcc2int (mM)) {
  V_T = r_T*(Mkcc2int)
  v_T = V_T*Kint*Cint/((Kint + R_T)*(Cint + R_T))
}

PROCEDURE driveforce(Cint (mM), Cout (mM), Kint (mM), Kout (mM)) {
  LOCAL T
  T = celsius + 273.15
  D = (1e3)*(R*T/F)*log((Kout + Cint)/(Kint + Cout))
  if (D > 0) {
    transport = v_T
  }
  if (D < 0) {
    transport = -v_T
  }
}
