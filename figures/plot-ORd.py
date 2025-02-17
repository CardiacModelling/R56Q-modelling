import matplotlib.pyplot as plt
import myokit
import numpy as np

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import argparse

parser = argparse.ArgumentParser(
    description='Plot ORd model results')
parser.add_argument("--RPR", action='store_true',
                    help="whether to use R56Q-RPR or not",
                    default=False)
parser.add_argument("--show", action='store_true',
                    help="whether to show figures instead of saving them",
                    default=False)
parser.add_argument("--ICaL", action='store_true',
                    help="whether to apply effect of RPR on ICaL",
                    default=False)
parser.add_argument('--dpi', type=int, default=100, metavar='N',
                    help='what DPI to use for figures (suggested 300 for publication quality)')
args = parser.parse_args()

RPR = args.RPR
ICaL = args.ICaL

WT = np.loadtxt('../cmaesfits/WT-model-M10-fit-staircase1-artefact.txt', unpack=True)
R56Q = np.loadtxt('../cmaesfits/R56Q-model-M10-fit-staircase1-artefact.txt', unpack=True)

# Get model
m = myokit.load_model('../model-and-protocols/Ord-original.mmt')

n_params = 12
parameters = []
for i in range(n_params):
    parameters.append('ikr.p'+str(i+1))

cell_mode = 0
fac = 1
R56Q_fac = 1
n_APs = 4

# Set cell type
m.set_value('cell.mode', cell_mode)
types = ['Endocardial', 'Epicardial', 'Midmyocardial']

# Block IKr
gkr = m.get('ikr.GKr').eval()
gcal = m.get('ical.PCa').eval()
m.get('ikr.GKr').set_rhs(fac*gkr)

# Create a pacing protocol
bcl = 400
p = myokit.pacing.blocktrain(bcl, 0.5, offset=20)

# Create simulation with APD tracking enabled
s = myokit.Simulation(m, p, apd_var='membrane.V')

vt = 0.9 * s.state()[m.get('membrane.V').indice()]

# Run, using the same threshold as before
s.pre(bcl * 100)
d, apds = s.run(bcl*n_APs, apd_threshold=vt)

m = myokit.parse_model('''
[[model]]
name: ohara-2011
date: 02/07/19
author: Dominic Whittaker
desc: """
    Myokit implementation of the O'Hara-Rudy dynamic (ORd) model for the
    undiseased human ventricular action potential and calcium transient.

    Based on the matlab code published on the rudylab website (see below) and
    checked against the supplement published with [1]

    References:
    [1] O'Hara et al. (2011) Simulation of the Undiseased Human Cardiac
    Ventricular Action Potential: Model Formulation and Experimental
    Validation. PLoS Computational Biology
    doi: 10.1371/journal.pcbi.1002061


    Original copyright notice:
    ---------------------------------------------------------------------------
    MATLAB Implementation of the O'Hara-Rudy dynamic (ORd) model for the
    undiseased human ventricular action potential and calcium transient

    The ORd model is described in the article "Simulation of the Undiseased
    Human Cardiac Ventricular Action Potential: Model Formulation and
    Experimental Validation"
    by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy

    The article and supplemental materails are freely available in the
    Open Access jounal PLoS Computational Biology
    Link to Article:
    http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061

    Email: tom.ohara@gmail.com / rudy@wustl.edu
    Web: http://rudylab.wustl.edu

    The ORd model is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. The ORd model is distributed in the hope that
    it will be useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details (www.gnu.org/licenses)
    """
# Initial values
membrane.V      = -87   # V=-87              1
sodium.Nai      = 7     # nai=7              2
sodium.Na_ss    = 7     # nass=nai           3
potassium.Ki    = 145   # ki=145             4
potassium.K_ss  = 145   # kss=ki             5
calcium.Cai     = 1e-4  # cai=1.0e-4         6
calcium.Ca_ss   = 1e-4  # cass=cai           7
calcium.Ca_nsr  = 1.2   # cansr=1.2          8
calcium.Ca_jsr  = 1.2   # cajsr=cansr        9
ina.m           = 0     # m=0               10
ina.hf          = 1     # hf=1              11
ina.hs          = 1     # hs=1              12
ina.j           = 1     # j=1               13
ina.hsp         = 1     # hsp=1             14
ina.jp          = 1     # jp=1              15
inal.m          = 0     # mL=0              16
inal.h          = 1     # hL=1              17
inal.hp         = 1     # hLp=1             18
ito.a           = 0     # a=0               19
ito.if          = 1     # iF=1              20
ito.is          = 1     # iS=1
ito.ap          = 0     # ap=0              22
ito.ifp         = 1     # iFp=1             23
ito.isp         = 1     # iSp=1
ical.d          = 0     # d=0               25
ical.ff         = 1     # ff=1              26
ical.fs         = 1     # fs=1
ical.fcaf       = 1     # fcaf=1            28
ical.fcas       = 1     # fcas=1
ical.jca        = 1     # jca=1             30
ical.nca        = 0     # nca=0
ical.ffp        = 1     # ffp=1             32
ical.fcafp      = 1     # fcafp=1           33
ikr.c1 = 0
ikr.c2 = 1
ikr.i = 0
ikr.ic1 = 0
ikr.ic2 = 0
ikr.o = 0
iks.x1          = 0     # xs1=0             36
iks.x2          = 0     # xs2=0
ik1.x           = 1     # xk1=1             38
ryr.Jrelnp      = 0     # Jrelnp=0
ryr.Jrelp       = 0     # Jrelp=0           40
camk.CaMKt      = 0     # CaMKt=0           41

#
# Engine variables
#
[engine]
time = 0 in [ms] bind time
pace = 0 bind pace

#
# Membrane potential
# Page 5
#
[membrane]
dot(V) = -(i_ion + stimulus.i_stim)
    label membrane_potential
    in [mV]
i_ion = (sodium.INa_tot
        + sodium.INa_ss_tot
        + calcium.ICa_tot
        + calcium.ICa_ss_tot
        + potassium.IK_tot
        + potassium.IK_ss_tot
        )
    label cellular_current
    in [uA/uF]

#
# Stimulus current
# Page 5
#
[stimulus]
i_stim = engine.pace * amplitude
amplitude = -80 [uA/uF]

#
# Cell geometry
# Page 6
#
[cell]
mode = 0
    desc: The type of cell. Endo = 0, Epi = 1, Mid = 2
L = 0.01 [cm] : Cell length
r = 0.0011 [cm] :  Cell radius
vcell = 1000 * 3.14 * r * r * L
    in [uL]
    desc: Cell volume
Ageo = 2*3.14 * r * r + 2 * 3.14 * r * L
    in [cm^2]
    desc: Geometric cell area
Acap = 2 * Ageo
    in [cm^2]
    desc: Capacitative membrane area
AF = Acap / phys.F
vmyo = 0.68 * vcell
    in [uL]
    desc: Volume of the cytosolic compartment
vnsr = 0.0552 * vcell
    in [uL]
    desc: Volume of the NSR compartment
vjsr = 0.0048 * vcell
    in [uL]
    desc: Volume of the JSR compartment
vss = 0.02 * vcell
    in [uL]
    desc: Volume of the Submembrane space near the T-tubules

#
# Physical constants
# Page 2
#
[phys]
R = 8314  [J/kmol/K] : Gas constant
T = 310   [K] : Temperature
F = 96485 [C/mol] : Faraday's constant
RTF  = R*T/F
FRT  = F/(R*T)
FFRT = F*F/(R*T)

#
# Extracellular concentrations
# Page 5
#
[extra]
Nao = 140 [mmol/L] : Extracellular Na+ concentration
Cao = 1.8 [mmol/L] : Extracellular Ca2+ concentration
Ko  = 5.4 [mmol/L] : Extracellular K+ concentration

#
# Reversal potentials
# Page 6
#
[nernst]
ENa = phys.RTF * log(extra.Nao / sodium.Nai)
    in [mV]
    desc: Reversal potential for Sodium currents
EK = phys.RTF * log(extra.Ko / potassium.Ki)
    in [mV]
    desc: Reversal potential for Potassium currents
PNaK = 0.01833
    desc: Permeability ratio K+ to Na+
EKs = phys.RTF * log((extra.Ko + PNaK * extra.Nao) / (potassium.Ki + PNaK * sodium.Nai))
    desc: Reversal potential for IKs
    in [mV]

#
# INa :: Fast Sodium current
# Page 6
#
# The fast sodium current is modelled using a Hodgkin-Huxley type formulation
# including activation (m), slow and fast components of inactivation (h) and
# recovery from inactivation (j). The slow component of inactivation and
# recovery from inactivation have an alternative formulation for CaMKII-
# phosphorylated channels.
#
[ina]
use membrane.V
tm  = 1 / (6.765 * exp((V + 11.64) / 34.77) + 8.552 * exp(-(V + 77.42) / 5.955))
    desc: Time constant for m-gate
    in [ms]
sm  = 1 / (1 + exp((V + 39.57) / -9.871))
    desc: Steady state value for m-gate
dot(m) = (sm - m) / tm
    desc: Activation gate for INa channels
sh  = 1 / (1 + exp((V + 82.90) / 6.086))
    desc: Steady-state value for h-gate
thf = 1 / (1.432e-5 * exp((V + 1.196) / -6.285) + 6.1490 * exp((V + 0.5096) / 20.27))
    desc: Time constant for fast development of inactivation in INa
    in [ms]
ths = 1 / (0.009794 * exp((V + 17.95) / -28.05) + 0.3343 * exp((V + 5.7300) / 56.66))
    desc: Time constant for slow development of inactivation in INa
    in [ms]
Ahf = 0.99 : Fraction of INa channels with fast inactivation
Ahs = 1.0 - Ahf : Fraction of INa channels with slow inactivation
dot(hf) = (sh - hf) / thf
    desc: Fast componennt of the inactivation gate for INa channels
dot(hs) = (sh - hs) / ths
    desc: Slow componennt of the inactivation gate for non-phosphorylated INa channels
h = Ahf * hf + Ahs * hs
    desc: Inactivation gate for INa
tj = 2.038 + 1 / (0.02136 * exp((V + 100.6) / -8.281) + 0.3052 * exp((V + 0.9941) / 38.45))
    desc: Time constant for j-gate in INa
    in [ms]
sj = sh
    desc: Steady-state value for j-gate in INa
dot(j) = (sj - j) / tj
    desc: Recovery from inactivation gate for non-phosphorylated INa channels
# Phosphorylated channels
thsp = 3 * ths
    desc: Time constant for h-gate of phosphorylated INa channels
    in [ms]
shsp = 1 / (1 + exp((V + 89.1) / 6.086))
    desc: Steady-state value for h-gate of phosphorylated INa channels
dot(hsp) = (shsp - hsp) / thsp
    desc: Slow componennt of the inactivation gate for phosphorylated INa channels
hp = Ahf * hf + Ahs * hsp
    desc: Inactivation gate for phosphorylated INa channels
tjp = 1.46 * tj
    desc: Time constant for the j-gate of phosphorylated INa channels
    in [ms]
dot(jp) = (sj - jp) / tjp
    desc: Recovery from inactivation gate for phosphorylated INa channels
# Current
GNa = 75 : Maximum conductance of INa channels
INa = GNa * (V - nernst.ENa) * m^3 * ((1 - camk.f) * h * j + camk.f * hp * jp)
    in [uA/uF]
    desc: Fast sodium current

#
# INaL :: Late component of the Sodium current
# Page 7
#
[inal]
use membrane.V
use ina.tm
sm = 1 / (1 + exp((V + 42.85) / -5.264))
    desc: Steady state value of m-gate for INaL
dot(m) = (sm - m) / tm
    desc: Activation gate for INaL
th = 200.0 [ms] : Time constant for inactivation of non-phosphorylated INaL channels
sh = 1 / (1 + exp((V + 87.61) / 7.488))
    desc: Steady-state value for inactivation of non-phosphorylated INaL channels
dot(h) = (sh - h) / th
    desc: Inactivation gate for non-phosphorylated INaL channels
thp = 3 * th
    in [ms]
    desc: Time constant for inactivation of phosphorylated INaL channels
shp = 1 / (1 + exp((V + 93.81) / 7.488))
    desc: Steady state value for inactivation of phosphorylated INaL channels
dot(hp) = (shp - hp) / thp
    desc: Inactivation gate for phosphorylated INaL channels
# Current
GNaL = 0.0075 : Maximum conductance of INaL
f_gnal = if(cell.mode == 1, 0.6, 1)
    desc: Adjustment for different cell types
INaL = f_gnal * GNaL * (V - nernst.ENa) * m * ((1 - camk.f) * h + camk.f * hp)

#
# Ito :: Transient outward Potassium current
# page 8
#
[ito]
use membrane.V
ta = 1.0515 / (one + two)
    one = 1 / (1.2089 * (1 + exp((V - 18.4099) / -29.3814)))
    two = 3.5 / (1 + exp((V + 100) / 29.3814))
    desc: Time constant for Ito activation
    in [ms]
sa = 1 / (1 + exp((V - 14.34) / -14.82))
    desc: Steady-state value for Ito activation
dot(a) = (sa - a) / ta
    desc: Ito activation gate
si = 1 / (1 + exp((V + 43.94) / 5.711))
    desc: Steady-state value for Ito inactivation
delta_epi = if(cell.mode == 1,
    1 - (0.95 / (1 + exp((V + 70) / 5))),
    1)
    desc: Adjustment for different cell types
tif = (4.562 + 1 / (0.3933 * exp((V+100) / -100) + 0.08004 * exp((V + 50) / 16.59))) * delta_epi
    desc: Time constant for fast component of Ito inactivation
    in [ms]
tis = (23.62 + 1 / (0.001416 * exp((V + 96.52) / -59.05)+ 1.780e-8 * exp((V + 114.1) / 8.079))) * delta_epi
    desc: Time constant for slow component of Ito inactivation
    in [ms]
dot(if) = (si - if) / tif
    desc: Fast component of Ito activation
dot(is) = (si - is) / tis
    desc: Slow component of Ito activation
Aif = 1 / (1 + exp((V - 213.6) / 151.2))
    desc: Fraction of fast inactivating Ito channels
Ais = 1 - Aif
    desc: Fraction of slow inactivating Ito channels
i = Aif * if + Ais * is
    desc: Inactivation gate for non-phosphorylated Ito
dot(ap) = (assp - ap) / ta
    assp=1.0/(1.0+exp((-(V-24.34))/14.82))
dti_develop = 1.354 + 1e-4 / (exp((V - 167.4) / 15.89) + exp((V - 12.23) / -0.2154))
dti_recover = 1 - 0.5 / (1 + exp((V+70) / 20))
tifp = dti_develop * dti_recover * tif
    desc: Time constant for fast component of inactivation of phosphorylated Ito channels
    in [ms]
tisp = dti_develop * dti_recover * tis
    desc: Time constant for slot component of inactivation of phosphorylated Ito channels
    in [ms]
dot(ifp) = (si - ifp) / tifp
    desc: Fast component of inactivation of phosphorylated Ito channels
dot(isp) = (si - isp) / tisp
    desc: Slow component of inactivation of phosphorylated Ito channels
ip = Aif * ifp + Ais * isp
    desc: Inactivation gate for phosphorylated Ito channels
# Current
Gto = if(cell.mode == 0, 0.02, 0.08)
    desc: Maximum conductance of Ito channels
Ito = Gto * (V - nernst.EK) * ((1 - camk.f) * a * i + camk.f * ap * ip)
    desc: Transient outward Potassium current

#
# ICaL  :: L-type Calcium current
# ICaNa :: Sodium current through the L-type Calcium channel
# ICaK  :: Potassium current through the L-type Calcium channel
# Page 9
#
# The ICaL channel is modeled using activation, inactivation (fast and slow),
# Ca-dependent inactivation (fast and slow) and recovery from Ca-dependent
# inactivation.
# Inactivation and Ca-dependent inactivation have an alternative formulation
# for CaMKII phosphorylated channels.
#
#
[ical]
use membrane.V
vf = V * phys.FRT
vff = V * phys.FFRT
# Activation
act_fac = 0
sd = 1 / (1 + exp((V + 3.94 + act_fac) / -4.23))
    desc: Steady-state value for activation gate of ICaL channel
td = 0.6 + 1 / (exp(-0.05 * (V + 6)) + exp(0.09 * (V + 14)))
    desc: Time constant for activation gate of ICaL channel
    in [ms]
dot(d) = (sd - d) / td
    desc: Activation gate of ICaL channel
# Inactivation
sf = 1 / (1 + exp((V + 19.58) / 3.696))
    desc: Steady-state value for inactivation gate of ICaL channel
tff = 7 + 1 / (0.0045 * exp((V + 20) / -10) + 0.0045 * exp((V + 20) / 10))
    desc: Time constant for fast inactivation of ICaL channels
    in [ms]
tfs = 1000 + 1 / (0.000035 * exp((V + 5) / -4) + 0.000035 * exp((V + 5) / 6))
    desc: Time constant for fast inactivation of ICaL channels
    in [ms]
dot(ff) = (sf - ff) / tff
    desc: Fast inactivation of ICaL channels
dot(fs) = (sf - fs) / tfs
    desc: Slow inactivation of ICaL channels
Aff = 0.6 : Fraction of ICaL channels with fast inactivation
Afs = 1 - Aff : Fraction of ICaL channels with slow inactivation
f = Aff * ff + Afs * fs
    desc: Inactivation of ICaL channels
# Ca-dependent inactivation
sfca = sf
    desc: Steady-state value for Ca-dependent inactivation of ICaL channels
tfcaf = 7 + 1 / (0.04 * exp((V - 4) / -7) + 0.04 * exp((V - 4) / 7))
    desc: Time constant for fast Ca-dependent inactivation of ICaL channels
    in [ms]
tfcas = 100 + 1 / (0.00012 * exp(V / -3) + 0.00012 * exp(V / 7))
    desc: Time constant for slow Ca-dependent inactivation of ICaL channels
    in [ms]
Afcaf = 0.3 + 0.6 / (1 + exp((V - 10) / 10))
    desc: Fraction of ICaL channels with fast Ca-dependent inactivation
Afcas = 1 - Afcaf
    desc: Fraction of ICaL channels with slow Ca-dependent inactivation
dot(fcaf) = (sfca - fcaf) / tfcaf
    desc: Fast Ca-dependent inactivation of ICaL channels
dot(fcas) = (sfca - fcas) / tfcas
    desc: Slow Ca-dependent inactivation of ICaL channels
fca = Afcaf * fcaf + Afcas * fcas
    desc: Ca-dependent inactivation of ICaL channels
# Recovery from Ca-dependent inactivation
tjca = 75 [ms] : Time constant of recovery from Ca-dependent inactivation
dot(jca) = (sfca - jca) / tjca
    desc: Recovery from Ca-dependent inactivation
# Inactivation of phosphorylated channels
tffp = 2.5 * tff
    in [ms]
    desc: Time constant for fast inactivation of phosphorylated ICaL channels
dot(ffp) = (sf - ffp) / tffp
    desc: Fast inactivation of phosphorylated ICaL channels
fp = Aff * ffp + Afs * fs
    desc: Inactivation of phosphorylated ICaL channels
# Ca-dependent inactivation of phosphorylated channels
tfcafp = 2.5 * tfcaf
    in [ms]
    desc: Time constant for fast Ca-dependent inactivation of phosphorylated ICaL channels
dot(fcafp) = (sfca - fcafp) / tfcafp
    desc: Fast Ca-dependent inactivation of phosphorylated ICaL channels
fcap = Afcaf * fcafp + Afcas * fcas
    desc: Ca-dependent inactivation of phosphorylated ICaL channels
# Fraction of channels in Ca-depdent inactivation mode
dot(nca) = anca * k2n - nca*km2n
    anca = 1 / (k2n / km2n + (1 + Kmn / calcium.Ca_ss)^4.0)
    Kmn = 0.002
    k2n = 1000
    km2n = jca * 1.0
    desc: Fraction of channels in Ca-depdent inactivation mode
# Total currents through ICaL channel
PhiCaL  = 4 * vff *(       calcium.Ca_ss  * exp(2 * vf) - 0.341 * extra.Cao) / (exp(2 * vf) - 1)
PhiCaNa = 1 * vff *(0.75 * sodium.Na_ss   * exp(1 * vf) - 0.75  * extra.Nao) / (exp(1 * vf) - 1)
PhiCaK  = 1 * vff *(0.75 * potassium.K_ss * exp(1 * vf) - 0.75  * extra.Ko ) / (exp(1 * vf) - 1)
PCa = piecewise(cell.mode == 0, base, cell.mode == 1, 1.2*base, 2.5*base)
    base = 0.0001
PCap   = 1.1      * PCa
PCaNa  = 0.00125  * PCa
PCaK   = 3.574e-4 * PCa
PCaNap = 0.00125  * PCap
PCaKp  = 3.574e-4 * PCap
g  = d * (f  * (1 - nca) + jca * fca  * nca)
    desc: Conductivity of non-phosphorylated ICaL channels
gp = d * (fp * (1 - nca) + jca * fcap * nca)
    desc: Conductivity of phosphorylated ICaL channels
ICaL   = (1 - camk.f) * PCa   * PhiCaL  * g + camk.f * PCap   * PhiCaL  * gp
    desc: L-type Calcium current
    in [uA/uF]
ICaNa  = (1 - camk.f) * PCaNa * PhiCaNa * g + camk.f * PCaNap * PhiCaNa * gp
    desc: Sodium current through ICaL channels
    in [uA/uF]
ICaK   = (1 - camk.f) * PCaK  * PhiCaK  * g + camk.f * PCaKp  * PhiCaK  * gp
    desc: Potassium current through ICaL channels
    in [uA/uF]

#
# IKr :: Rapid delayed rectifier Potassium current
# Page 11
[ikr]
use membrane.V

dot(c1) = a1 * c2 + ah * ic1 + b2 * o - (b1 + bh + a2) * c1
dot(c2) = b1 * c1 + ah * ic2 - (a1 + bh) * c2
dot(i) = a2 * ic1 + bh * o - (b2 + ah) * i
dot(ic1) = a1 * ic2 + bh * c1 + b2 * i - (b1 + ah + a2) * ic1
dot(ic2) = b1 * ic1 + bh * c2 - (ah + a1) * ic2
dot(o) = a2 * c1 + ah * i - (b2 + bh) * o

a1 = p1 * exp(p2 * V)
b1 = p3 * exp(-p4 * V)
bh = p5 * exp(p6 * V)
ah = p7 * exp(-p8 * V)
a2 = p9 * exp(p10 * V)
b2 = p11 * exp(-p12 * V)

p1 = 2.26e-4 [1/ms]
p2 = 0.06990 [1/mV]
p3 = 3.45e-5 [1/ms]
p4 = 0.05462 [1/mV]
p5 = 0.08730 [1/ms]
p6 = 8.91e-3 [1/mV]
p7 = 5.15e-3 [1/ms]
p8 = 0.03158 [1/mV]
p9 = 0.08730 [1/ms]
p10 = 8.91e-3 [1/mV]
p11 = 5.15e-3 [1/ms]
p12 = 0.03158 [1/mV]
p13 = 0.05 [mS/uF]

base = p13
GKr = piecewise(cell.mode == 0, base, cell.mode == 1, 1.3 * base, 0.8 * base)
    #base = p13

IKr = GKr * o * (V - nernst.EK)
    desc: Rapid delayed Potassium current
    in [uA/uF]

# IKs :: Slow delayed rectifier Potassium current
# Page 11
#
# Modelled with two activation channels
#
[iks]
use membrane.V
sx  = 1 / (1 + exp((V + 11.60) / -8.932))
    desc: Steady-state value for activation of IKs channels
tx1 = 817.3 + 1 / (2.326e-4 * exp((V + 48.28) / 17.80) + 0.001292 * exp((V + 210) / -230))
    desc: Time constant for slow, low voltage IKs activation
dot(x1) = (sx - x1) / tx1
    desc: Slow, low voltage IKs activation
tx2 = 1 / (0.01 * exp((V - 50) / 20) + 0.0193 * exp((V + 66.54) / -31))
    desc: Time constant for fast, high voltage IKs activation
dot(x2) = (sx - x2) / tx2
    desc: Fast, high voltage IKs activation
KsCa = 1 + 0.6 / (1.0 + (3.8e-5 / calcium.Cai)^1.4)
    desc: Maximum conductance for IKs
GKs = if(cell.mode == 1, 1.4 * base, base)
    base = 0.0034
    desc: Conductivity adjustment for cell type
IKs = GKs * KsCa * x1 * x2 * (V - nernst.EKs)
    desc: Slow delayed rectifier Potassium current

#
# IK1 :: Inward rectifier Potassium current
# Page 12
#
# Modelled with an activation channel and an instantaneous inactivation channel
#
[ik1]
use membrane.V
sx = 1 / (1 + exp(-(V + 2.5538 * extra.Ko + 144.59) / (1.5692 * extra.Ko + 3.8115)))
    desc: Steady-state value for activation of IK1 channels
tx = 122.2 / (exp((V + 127.2) / -20.36) + exp((V + 236.8) / 69.33))
    desc: Time constant for activation of IK1 channels
dot(x) = (sx - x) / tx
    desc: Activation of IK1 channels
r = 1 / (1 + exp((V + 105.8 - 2.6 * extra.Ko) / 9.493))
    desc: Inactivation of IK1 channels
GK1 = piecewise(cell.mode == 0, base, cell.mode == 1, 1.2*base, 1.3*base)
    base = 0.1908
    desc: Conductivity of IK1 channels, cell-type dependent
IK1 = GK1 * sqrt(extra.Ko) * r * x * (V - nernst.EK)
    desc: Inward rectifier Potassium current

#
# INaCa :: Sodium/Calcium exchange current
# page 12
#
[inaca]
use membrane.V
use extra.Nao, extra.Cao
use sodium.Nai, calcium.Cai
kna1   = 15.0
kna2   = 5.0
kna3   = 88.12
kasymm = 12.5
wna    = 6.0e4
wca    = 6.0e4
wnaca  = 5.0e3
kcaon  = 1.5e6
kcaoff = 5.0e3
qna    = 0.5224
qca    = 0.1670
hca    = exp(qca * V * phys.FRT)
hna    = exp(qna * V * phys.FRT)
# Parameters h
h1  = 1 + Nai / kna3 * (1 + hna)
h2  = (Nai * hna) / (kna3 * h1)
h3  = 1 / h1
h4  = 1 + Nai / kna1 * (1 + Nai / kna2)
h5  = Nai * Nai / (h4 * kna1 * kna2)
h6  = 1 / h4
h7  = 1 + Nao / kna3 * (1 + 1 / hna)
h8  = Nao / (kna3 * hna * h7)
h9  = 1 / h7
h10 = kasymm + 1 + Nao / kna1 * (1 + Nao / kna2)
h11 = Nao * Nao / (h10 * kna1 * kna2)
h12 = 1 / h10
# Parameters k
k1   = h12 * Cao * kcaon
k2   = kcaoff
k3p  = h9 * wca
k3pp = h8 * wnaca
k3   = k3p + k3pp
k4p  = h3 * wca / hca
k4pp = h2 * wnaca
k4   = k4p + k4pp
k5   = kcaoff
k6   = h6 * Cai * kcaon
k7   = h5 * h2 * wna
k8   = h8 * h11 * wna
x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
E1 = x1 / (x1 + x2 + x3 + x4)
E2 = x2 / (x1 + x2 + x3 + x4)
E3 = x3 / (x1 + x2 + x3 + x4)
E4 = x4 / (x1 + x2 + x3 + x4)
KmCaAct = 150.0e-6
allo    = 1 / (1 + (KmCaAct / Cai)^2.0)
JncxNa  = 3 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
JncxCa  = E2 * k2 - E1 * k1
Gncx = piecewise(cell.mode == 0, base, cell.mode == 1, 1.1*base, 1.4*base)
    base = 0.0008
INaCa = 0.8 * Gncx * allo * (JncxNa + 2 * JncxCa)
    desc: Sodium/Calcium exchange current
    in [uA/uF]

#
# INaCa_ss :: Sodium/Calcium exchanger current into the L-type subspace
# Page 12
#
[inacass]
use membrane.V
use extra.Nao, extra.Cao
use sodium.Na_ss, calcium.Ca_ss
h1  = 1 + Na_ss / inaca.kna3 * (1 + inaca.hna)
h2  = (Na_ss * inaca.hna)/(inaca.kna3 * h1)
h3  = 1 / h1
h4  = 1 + Na_ss / inaca.kna1 * (1 + Na_ss / inaca.kna2)
h5  = Na_ss * Na_ss /(h4 * inaca.kna1 * inaca.kna2)
h6  = 1 / h4
h7  = 1 + Nao / inaca.kna3 * (1 + 1 / inaca.hna)
h8  = Nao / (inaca.kna3 * inaca.hna * h7)
h9  = 1/h7
h10 = inaca.kasymm + 1 + Nao / inaca.kna1 * (1 + Nao / inaca.kna2)
h11 = Nao * Nao / (h10 * inaca.kna1 * inaca.kna2)
h12 = 1/h10
k1   = h12 * Cao * inaca.kcaon
k2   = inaca.kcaoff
k3p  = h9 * inaca.wca
k3pp = h8 * inaca.wnaca
k3   = k3p + k3pp
k4p  = h3 * inaca.wca / inaca.hca
k4pp = h2 * inaca.wnaca
k4   = k4p + k4pp
k5   = inaca.kcaoff
k6   = h6 * Ca_ss * inaca.kcaon
k7   = h5 * h2 * inaca.wna
k8   = h8 * h11 * inaca.wna
x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3)
x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8)
x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3)
x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8)
E1 = x1 / (x1 + x2 + x3 + x4)
E2 = x2 / (x1 + x2 + x3 + x4)
E3 = x3 / (x1 + x2 + x3 + x4)
E4 = x4 / (x1 + x2 + x3 + x4)
KmCaAct = 150.0e-6
allo    = 1 / (1 + (KmCaAct / Ca_ss)^2)
JncxNa  = 3 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp
JncxCa  = E2 * k2 - E1 * k1
INaCa_ss = 0.2 * inaca.Gncx * allo * (JncxNa + 2 * JncxCa)
    desc: Sodium/Calcium exchange current into the T-Tubule subspace
    in [uA/uF]

#
# INaK :: Sodium/Potassium ATPase current
# Page 14
#
[inak]
use membrane.V
use extra.Nao, sodium.Nai, sodium.Na_ss
use extra.Ko, potassium.Ki, potassium.K_ss
k1p = 949.5
k1m = 182.4
k2p = 687.2
k2m = 39.4
k3p = 1899.0
k3m = 79300.0
k4p = 639.0
k4m = 40.0
Knai0 = 9.073
Knao0 = 27.78
delta = -0.1550
Knai = Knai0 * exp(delta * V * phys.FRT / 3)
Knao = Knao0 * exp((1.0-delta) * V * phys.FRT / 3)
Kki    = 0.5
Kko    = 0.3582
MgADP  = 0.05
MgATP  = 9.8
Kmgatp = 1.698e-7
H      = 1.0e-7
eP     = 4.2
Khp    = 1.698e-7
Knap   = 224.0
Kxkur  = 292.0
P = eP / (1 + H / Khp + Nai / Knap + Ki / Kxkur)
a1 = (k1p * (Nai / Knai)^3) / ((1 + Nai / Knai)^3 + (1 + Ki / Kki)^2 - 1)
b1 = k1m * MgADP
a2 = k2p
b2 = (k2m * (Nao / Knao)^3) / ((1 + Nao / Knao)^3 + (1 + Ko / Kko)^2 - 1)
a3 = (k3p * (Ko / Kko)^2  ) / ((1 + Nao / Knao)^3 + (1 + Ko / Kko)^2 - 1)
b3 = (k3m * P * H)/(1 + MgATP / Kmgatp)
a4 = (k4p * MgATP / Kmgatp) / (1 + MgATP / Kmgatp)
b4 = (k4m * (Ki / Kki)^2) / ((1 + Nai / Knai)^3 + (1 + Ki / Kki)^2 - 1)
x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2
x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4
x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1
x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1
E1 = x1 / (x1 + x2 + x3 + x4)
E2 = x2 / (x1 + x2 + x3 + x4)
E3 = x3 / (x1 + x2 + x3 + x4)
E4 = x4 / (x1 + x2 + x3 + x4)
JnakNa = 3 * (E1 * a3 - E2 * b3)
JnakK  = 2 * (E4 * b1 - E3 * a1)
Pnak = piecewise(cell.mode == 0, base, cell.mode == 1, 0.9*base, 0.7*base)
    base = 30
INaK = Pnak * (JnakNa + JnakK)
    desc: Sodium/Potassium ATPase current
    in [uA/uF]

#
# IKb :: Background Potassium current
# Page 15
#
[ikb]
use membrane.V
xkb = 1 / (1 + exp((V - 14.48) / -18.34))
GKb = if(cell.mode == 1, 0.0018, 0.003)
IKb = GKb * xkb * (V - nernst.EK)
    desc: Background Potassium current
    in [uA/uF]

#
# INab :: Background Sodium current
# Page 15
#
[inab]
use membrane.V
PNab = 3.75e-10
INab = PNab * V * phys.FFRT * (sodium.Nai * evf - extra.Nao) / (evf - 1)
    evf = exp(V * phys.FRT)
    desc: Background Sodium current
    in [uA/uF]

#
# ICab :: Background Calcium current
# Page 15
#
[icab]
use membrane.V
PCab=2.5e-8
ICab = PCab * 4 * V * phys.FFRT * (calcium.Cai * evf2 - 0.341 * extra.Cao) / (evf2 - 1)
    evf2 = exp(2 * V * phys.FRT)
    desc: Background Calcium current
    in [uA/uF]

#
# IpCa :: Sarcolemmal Calcium pump current
# Page 15
#
[ipca]
use membrane.V
GpCa = 0.0005
IpCa = GpCa * calcium.Cai / (0.0005 + calcium.Cai)
    desc: Sarcolemmal Calcium pump current
    in [uA/uF]

#
# Jrel :: SR Calcium release flux via Ryanodine receptor
# Page 17
#
[ryr]
use membrane.V
bt=4.75
a_rel=0.5*bt
Jrel_inf = if(cell.mode == 2, 1.7 * base, base)
    base = a_rel * -ical.ICaL / (1 + (1.5 / calcium.Ca_jsr)^8)
dot(Jrelnp) = (Jrel_inf - Jrelnp) / tau_rel
    tau_rel = if(value < 0.001, 0.001, value)
    value = bt / (1.0 + 0.0123 / calcium.Ca_jsr)
btp = 1.25*bt
a_relp = 0.5*btp
Jrel_infp = if(cell.mode == 2, 1.7*base, base)
    base = a_relp * -ical.ICaL / (1 + (1.5 / calcium.Ca_jsr)^8)
dot(Jrelp) = (Jrel_infp - Jrelp) / tau_relp
    tau_relp = if(value < 0.001, 0.001, value)
    value = btp / (1 + 0.0123 / calcium.Ca_jsr)
Jrel = (1 - camk.f) * Jrelnp + camk.f * Jrelp
    desc: SR Calcium release flux via Ryanodine receptor
    in [mmol/L/ms]

#
# Jup :: Calcium uptake via SERCA pump
# Page 17
#
[serca]
use calcium.Cai, calcium.Ca_jsr, calcium.Ca_nsr
f = if(cell.mode == 1, 1.3, 1)
Jupnp = f * (0.004375 * Cai / (Cai + 0.00092))
Jupp  = f * (2.75 * 0.004375 * Cai / (Cai + 0.00092 - 0.00017))
Jleak = 0.0039375 * Ca_nsr / 15
    in [mmol/L/ms]
Jup = (1 - camk.f) * Jupnp + camk.f * Jupp - Jleak
    desc: Total Ca2+ uptake, via SERCA pump, from myoplasm to nsr
    in [mmol/L/ms]
Jtr = (Ca_nsr - Ca_jsr) / 100
    desc: Ca2+ translocation from nsr to jsr
    in [mmol/L/ms]

#
# Diffusion fluxes
# Page 16
#
[diff]
JdiffNa = (sodium.Na_ss - sodium.Nai) / 2
JdiffK  = (potassium.K_ss  - potassium.Ki)  / 2
Jdiff   = (calcium.Ca_ss - calcium.Cai) / 0.2

#
# Intracellular Sodium concentrations
# Page 18
#
[sodium]
use cell.AF, cell.vss, cell.vmyo
INa_tot    = ina.INa + inal.INaL + inab.INab + 3*inaca.INaCa + 3*inak.INaK
dot(Nai)   = -INa_tot * AF / vmyo + diff.JdiffNa * vss / vmyo
    desc: Intracellular Potassium concentration
INa_ss_tot = ical.ICaNa + 3*inacass.INaCa_ss
dot(Na_ss) = -INa_ss_tot * AF / vss - diff.JdiffNa

#
# Intracellular Potassium concentrations
# Page 18
#
[potassium]
use cell.AF, cell.vss, cell.vmyo
IK_tot = (ito.Ito
        + ikr.IKr
        + iks.IKs
        + ik1.IK1
        + ikb.IKb
        - 2 * inak.INaK
        )
IK_ss_tot = ical.ICaK
dot(Ki)  = -(IK_tot + stimulus.i_stim) * AF / vmyo + diff.JdiffK * vss / vmyo
    desc: Intracellular Potassium concentration
dot(K_ss) = -IK_ss_tot * AF / vss - diff.JdiffK
    desc: Potassium concentration in the T-Tubule subspace

#
# Intracellular Calcium concentrations and buffers
# Page 18
#
[calcium]
use cell.AF, cell.vmyo, cell.vnsr, cell.vjsr, cell.vss
cmdnmax = if(cell.mode == 1, 1.3*base, base)
    base = 0.05
kmcmdn  = 0.00238
trpnmax = 0.07
kmtrpn  = 0.0005
BSRmax  = 0.047
KmBSR   = 0.00087
BSLmax  = 1.124
KmBSL   = 0.0087
csqnmax = 10.0
kmcsqn  = 0.8
ICa_tot = ipca.IpCa + icab.ICab - 2*inaca.INaCa
dot(Cai)  = buff * (-ICa_tot * AF / (2*vmyo) - serca.Jup  * vnsr / vmyo + diff.Jdiff * vss / vmyo )
    buff = 1 / (1 + cmdnmax * kmcmdn / (a*a) + trpnmax * kmtrpn / (b*b))
    a = kmcmdn + calcium.Cai
    b = kmtrpn + calcium.Cai
    desc: Intracellular Calcium concentratium
    in [mmol/L]
ICa_ss_tot = ical.ICaL - 2 * inacass.INaCa_ss
dot(Ca_ss) = buff * (-ICa_ss_tot * AF / (2*vss) + ryr.Jrel * vjsr / vss - diff.Jdiff )
    buff = 1 / (1 + BSRmax * KmBSR / (a*a) + BSLmax * KmBSL / (b*b))
    a = KmBSR + calcium.Ca_ss
    b = KmBSL + calcium.Ca_ss
    desc: Calcium concentratium in the T-Tubule subspace
    in [mmol/L]
dot(Ca_jsr) = buff * (serca.Jtr - ryr.Jrel)
    buff = 1 / (1 + csqnmax * kmcsqn / (a * a))
    a = kmcsqn + calcium.Ca_jsr
    desc: Calcium concentration in the JSR subspace
    in [mmol/L]
dot(Ca_nsr) = serca.Jup - serca.Jtr * vjsr / vnsr
    desc: Calcium concentration in the NSR subspace
    in [mmol/L]

#
# CaMKII signalling
#
[camk]
KmCaMK = 0.15
aCaMK  = 0.05
bCaMK  = 0.00068
CaMKo  = 0.05
KmCaM  = 0.0015
CaMKb  = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / calcium.Ca_ss)
CaMKa  = CaMKb + CaMKt
dot(CaMKt) = aCaMK * CaMKb * CaMKa - bCaMK * CaMKt
f = 1 / (1 + KmCaMK / CaMKa)
    desc: Fraction of phosphorylated channels'''
    )
m.set_value('cell.mode', cell_mode)
types = ['Endocardial', 'Epicardial', 'Midmyocardial']

# WT
m.get('ikr.base').set_rhs(fac*0.038) 
s = myokit.Simulation(m, p, apd_var='membrane.V')

for i, name in enumerate(parameters):
    if not RPR:
        s.set_constant(name, WT[i])
    else:
        s.set_constant(name, R56Q[i])

vt = 0.9 * s.state()[m.get('membrane.V').indice()]

# Run, using the same threshold as before
s.pre(bcl * 100)
d2, apds = s.run(bcl*n_APs, apd_threshold=vt)

# R56Q
m.get('ikr.base').set_rhs(fac*R56Q_fac*0.038)
if RPR and ICaL:
    m.get('ical.PCa').set_rhs(gcal*0.9)
s = myokit.Simulation(m, p, apd_var='membrane.V')

for i, name in enumerate(parameters):
    if not RPR:
        s.set_constant(name, R56Q[i])
    else:
        s.set_constant(name, R56Q[i+n_params])

vt = 0.9 * s.state()[m.get('membrane.V').indice()]

# Run, using the same threshold as before
s.pre(bcl * 100)
d3, apds = s.run(bcl*n_APs, apd_threshold=vt)

if not RPR:
    colors = [(0.11765,0.23529,1.0), 'red']
else:
    colors = ['red', 'orange']

# Display the result
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 2.5), dpi=args.dpi, constrained_layout=True)

ax1.set_ylabel('Voltage (mV)', fontsize=9)
ax1.plot(d2['engine.time'], d2['membrane.V'], linewidth=1, c=colors[0])
ax1.plot(d3['engine.time'], d3['membrane.V'], linewidth=1, c=colors[1])
[label.set_visible(False) for label in ax1.get_xticklabels()]
if not RPR:
    ax1.legend(['WT', 'R56Q'], fontsize=8, ncol=2, loc='upper right')
else:
    ax1.legend(['R56Q', 'R56Q-RPR'], fontsize=8, ncol=2, loc='upper right')
ax1.grid(True)

ax2.plot(d2['engine.time'], d2['ikr.IKr'], linewidth=1, c=colors[0])
ax2.plot(d3['engine.time'], d3['ikr.IKr'], linewidth=1, c=colors[1])
ax2.set_ylabel('Current (A/F)', fontsize=9)
ax2.set_xlabel('Time (ms)', fontsize=9)
ax2.grid(True)

axins = inset_axes(ax2, 1,1, loc='center right')
axins.grid(True)
[label.set_visible(False) for label in axins.get_xticklabels()]
[label.set_visible(False) for label in axins.get_yticklabels()]
axins.plot(d2['engine.time'], d2['ikr.IKr'], linewidth=1, c=colors[0])
axins.plot(d3['engine.time'], d3['ikr.IKr'], linewidth=1, c=colors[1])

if not RPR:
    x1, x2, y1, y2 = 800, 850, -0.05, 0.85
else:
    x1, x2, y1, y2 = 800, 850, -0.05, 0.75

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

mark_inset(ax2, axins, loc1=2, loc2=3, fc="none", ec="0.5")

if args.show == True:
    plt.show()
else:
    filename = 'ORd-model-RPR-' + str(RPR) + '-ICaL-' + str(ICaL)
    plt.savefig('PNG_figures/' + filename + '.png')
