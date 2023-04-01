import ROOT
import os, sys
from math import *
import ctypes

path = os.getenv('GARFIELD_INSTALL')
if sys.platform == 'darwin':
  ROOT.gSystem.Load(path + '/lib/libmagboltz.dylib')
  ROOT.gSystem.Load(path + '/lib/libGarfield.dylib')
else:
  ROOT.gSystem.Load(path + '/lib/libmagboltz.so')
  ROOT.gSystem.Load(path + '/lib/libGarfield.so')

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('ar_80_co2_20_0T.gas')

efields = ROOT.std.vector('double')()
bfields = ROOT.std.vector('double')()
angles = ROOT.std.vector('double')()
gas.GetFieldGrid(efields, bfields, angles)
nE = efields.size()
print('  E [V/cm]  vE [cm/us]  alpha [1/cm] eta [1/cm]')
for i in range(nE):
  ve = ctypes.c_double(0.)
  gas.GetElectronVelocityE(i, 0, 0, ve)
  # Convert from cm/ns to cm/us.
  ve = 1.e3 * ve.value
  alpha = ctypes.c_double(0.)
  gas.GetElectronTownsend(i, 0, 0, alpha)
  alpha = exp(alpha.value)
  eta = ctypes.c_double(0.)
  gas.GetElectronAttachment(i, 0, 0, eta)
  eta = exp(eta.value)
  print('{0:10.2f} {1:10.2f} {2:10.2f} {3:10.2f}'.format(efields[i], ve, alpha, eta))
