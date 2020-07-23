#!/usr/bin/python3
import sys
import os
import numpy as np
from Motor import *

if ((len(sys.argv) ==2) and ( (sys.argv[1]=='-h') or (sys.argv[1]=='-help'))):
  print(sys.argv[0]," -h,-help           : display help message")
  print(sys.argv[0],"  ")
  print("      -f : Current frequency")
  print("      -nu : Main current frequency")
  print("      -nu0 : initial motor speed")
  print("      -Gamma : Load friction")
  print("      -Gamma3 : Load friction cubic coefficient")
  print("      -I : Moment of inertia")
  print("      -R : Resistance of loop")
  print("      -Phi : loop maximum flux (A B)")
  print("      -Kb : back EM constant: : A*k/(L*R)")
  print("      -N : Number of loops")
  print("      -Ns : number of stator coils (6)")
  print("      -trel : Relaxation time")
  print("      -tmax : max time")
  print("      -sf : sampling frequency")
  print("      -esp0 : anisotropy factor for loop 0")
  print("      -esp_pre : rotor presession amp")
  print("      -nu_pre : rotor presession frequency")
  print("      -esp_I : Inertia presession amp")
  print("      -nu_I : Inertia presession frequency")
  print("      -esp_G : Gamma presession amp")
  print("      -nu_G : Gamma presession frequency")
  print("      -wfG : wave form for Gamma : 0:sin, 1:sin^3, 2 : square")
  sys.exit(0)

def_nu = nu = 50
def_nu0 = nu0 = nu*0.95
def_Gamma = Gamma = 0.2
def_Gamma3 = Gamma3 = 0
def_R = R = 1
def_I = I = 1
def_Phi = Phi = 1
def_Kb = Kb = 1
def_N = N = 16
def_Ns = Ns = 6
def_trel = trel = 5
def_tmax = tmax = 15
def_sf = sf = 1000
def_eps0 = eps0 = 0.0
def_eps_pre = eps_pre = 0.0
def_nu_pre = nu_pre = 0.0
def_eps_I = eps_I = 0.0
def_nu_I = nu_I = 0.0
def_eps_G = eps_G = 0.0
def_nu_G = nu_G = 0.0
def_wfG = wfG = 0.0

# Change the value here not above
#nu = 50
#nu0 = nu*0.95
#Gamma = 0.2
#Gamma3 = 0
#R = 1
#I = 1
#Phi = 1
#Kb = 1
#N = 16
#Ns = 6
#trel = 5
#tmax = 15
#sf = 1000
#eps0 = 0.0
#eps_pre = 0.0
#nu_pre = 0.0
#eps_I = 0.0
#nu_I = 0.0
#eps_G = 0.0
#nu_G = 0.0
#wfG=0

argi = 1
while(argi < len(sys.argv)):
  w = sys.argv[argi].split("=")
  if(w[0] == "-nu"):
    nu = float(w[1])
  elif(w[0] == "-nu0"):
    nu0 = float(w[1])
  elif(w[0] == "-Gamma"):
    Gamma = float(w[1])
  elif(w[0] == "-Gamma3"):
    Gamma3 = float(w[1])
  elif(w[0] == "-I"):
    I = float(w[1])
  elif(w[0] == "-Phi"):
    Phi = float(w[1])
  elif(w[0] == "-Kb"):
    Kb = float(w[1])
  elif(w[0] == "-N"):
    N = int(w[1])
  elif(w[0] == "-Ns"):
    Ns = int(w[1])
  elif(w[0] == "-trel"):
    trel = float(w[1])
  elif(w[0] == "-tmax"):
    tmax = float(w[1])
  elif(w[0] == "-sf"):
    sf = float(w[1])
  elif(w[0] == "-eps0"):
    eps0 = float(w[1])
  elif(w[0] == "-eps_pre"):
    eps_pre = float(w[1])
  elif(w[0] == "-nu_pre"):
    nu_pre = float(w[1])
  elif(w[0] == "-eps_I"):
    eps_I = float(w[1])
  elif(w[0] == "-nu_I"):
    nu_I = float(w[1])
  elif(w[0] == "-eps_G"):
    eps_G = float(w[1])
  elif(w[0] == "-nu_G"):
    nu_G = float(w[1])
  elif(w[0] == "-wfG"):
    wfG = int(w[1])
  else:
    print("Invalid parameter ", sys.argv[argi])
  argi +=1
Omega = nu*2*np.pi

# Make filename using only parameters with value different from the default.  
def mk_filename():
  f_name = "motor"
  if(abs(def_nu - nu) > 1e-10):
    f_name +="_nu{}".format(nu)
  if(abs(def_nu0 - nu0) > 1e-10):
    f_name +="_nu0{}".format(nu0)
  if(abs(def_Gamma - Gamma) > 1e-10):
    f_name +="_G{}".format(Gamma)
  if(abs(def_Gamma3 - Gamma3) > 1e-10):
    f_name +="_Gt{}".format(Gamma3)
  if(abs(def_I - I) > 1e-10):
    f_name +="_I{}".format(I)
  if(abs(def_R - R) > 1e-10):
    f_name +="_R{}".format(R)
  if(abs(def_Phi - Phi) > 1e-10):
    f_name +="_Phi{}".format(Phi)
  if(abs(def_Kb - Kb) > 1e-10):
    f_name +="_Kb{}".format(Kb)
  if(abs(def_N - N) > 1e-10):
    f_name +="_N{}".format(N)
  if(abs(def_Ns - Ns) > 1e-10):
    f_name +="_Ns{}".format(Ns)
  if(abs(def_trel - trel) > 1e-10):
    f_name +="_trel{}".format(trel)
  if(abs(def_tmax - tmax) > 1e-10):
    f_name +="_tmax{}".format(tmax)
  if(abs(def_sf - sf) > 1e-10):
    f_name +="_sf{}".format(sf)
  if(abs(def_eps0 - eps0) > 1e-10):
    f_name +="_eps0{}".format(eps0)
  if(abs(def_eps_pre - eps_pre) > 1e-10):
    f_name +="_eps_pre{}".format(eps_pre)
  if(abs(def_nu_pre - nu_pre) > 1e-10):
    f_name +="_nu_pre{}".format(nu_pre)
  if(abs(def_eps_I - eps_I) > 1e-10):
    f_name +="_eps_I{}".format(eps_I)
  if(abs(def_nu_I - nu_I) > 1e-10):
    f_name +="_nu_I{}".format(nu_I)
  if(abs(def_eps_G - eps_G) > 1e-10):
    f_name +="_eps_G{}".format(eps_G)
  if(abs(def_nu_G - nu_G) > 1e-10):
    f_name +="_nu_G{}".format(nu_G)
  if(abs(def_wfG - wfG) > 1e-10):
    f_name +="_wfG{}".format(wfG)
  f_name += ".txt"
  return(f_name)

fname = mk_filename()
print(fname)

#fname = "motor_O{}_G{}_I{}_Kt{}_Kb{}_N{}_Ns{}_trel{}_tmax{}_sf{}_eps0{}_eps_pre{}_nu_pre{}.txt".format(Omega, Gamma, I, Kt, Kb, N, Ns, trel, tmax, sf, eps0, eps_pre, nu_pre)
  
mot = Motor(V0=[0,nu0*2*np.pi], dt=1/sf , t0=0)
mot.set(nu, Gamma, Gamma3 ,I, R,  Phi, Kb, N, Ns, trel, tmax, sf, eps0, eps_pre, nu_pre, eps_I, nu_I, eps_G, nu_G, wfG) 
mot.iterate(tmax, fig_dt=-1)
mot.save(fname, trel)
print("omega/Omega=", mot.V_list[-1][1]/mot.Omega)
print("slip=", (mot.Omega-mot.V_list[-1][1])/(2*np.pi),"Hz")
