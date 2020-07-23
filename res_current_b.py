#!/usr/bin/python3
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Condor queue
if ((len(sys.argv) ==2) and ( (sys.argv[1]=='-h') or (sys.argv[1]=='-help'))):
  print(sys.argv[0]," -h,-help           : display help message")
  print(sys.argv[0],"  ")
  print("      -f : datafile")
  print("      -fr : reference spectral datafile")
  print("      -t0 : relaxation time")
  print("      -p : No of point for FFT = 2^p")
  print("      -nu_max : Max Frequency for FFT plot")
  sys.exit(0)

fname = ""
fname_ref = ""
t0 = 0
#p=13 # 8192 points
#p=15 # 32768 points
p=18 # 262144 points
nu_max = -1000

argi = 1
while(argi < len(sys.argv)):
  w = sys.argv[argi].split("=")
  if(w[0] == "-f"):
    fname = w[1]
  elif(w[0] == "-fr"):
    fname_ref = w[1]
  elif(w[0] == "-t0"):
    t0 = float(w[1])
  elif(w[0] == "-p"):
    p = int(w[1])
  elif(w[0] == "-nu_max"):
    nu_max = int(w[1])
  else:
    print("Invalid parameter ", sys.argv[argi])
  argi +=1

def read_data(fname):
  t_list=[]
  v_list=[]
  n = 0
  omega=0
  for line in open(fname):
    s = line.split(" ")
    t=float(s[0])
    if(t>= t0) and (n < N):
      t_list.append(t)
      v_list.append(float(s[3])) # 4th column is residual current
      omega += float(s[2])
      n +=1
      
  return(t_list, v_list, 50-omega/(n*2*np.pi))

def fft_data(fname, N):
  t_list, v_list,slip = read_data(fname) 

  print(fname," slip=", slip)
  tmax = t_list[-1]-t_list[0]
  #print(t_list[-1],t_list[0])
  freq = np.fft.fftfreq(N,tmax/N)
  freq = np.abs(freq)[0:N//2+1]

  print("N=",N)
  ham = np.hamming(N)
  #freq, psd = signal.welch(v_list[:N]*ham)
  
  v = np.fft.fft(v_list[:N]*ham) # FFT of f
  #psd = np.abs(v[0:N//2+1])*2/N/freq
  modfft = np.abs(v[0:N//2+1])*2/N
  #modfft[0] /= 2     # constant term
  #modfft[N//2] /= 2  # higher frequency
  return(freq, modfft, tmax)

  # No of point for FFT as power of 2
N=1
for i in range(p):
 N *=2 

freq, modfft, tmax = fft_data(fname, N)

"""t_list, v_list = read_data(fname) 
tmax = t_list[-1]-t_list[0]
freq = np.fft.fftfreq(N,tmax/N)
freq = np.abs(freq)[0:N//2+1]
v = np.fft.fft(v_list) # FFT of f
modfft = np.abs(v[0:N//2+1])*2/N
"""

nfft  = N
if(nu_max > 0):
  fmax = N/tmax
  nfft = int(N*nu_max/fmax)
  print(nfft,N,nu_max,tmax,fmax)
if(nfft> N): nfft = N

if(fname_ref != ""):
  freq_r, modfft_r, tmax_r = fft_data(fname_ref, N)
  val = modfft[:nfft]-modfft_r[:nfft]
  plt.xlabel(r'$\nu$')
  plt.ylabel(r'$|I|-|I_{ref}|$')
  plt.margins(0.1, 0.1)
  plt.title("Residual current spectrum")
  plt.plot(freq[:nfft],val,"b.") 
  plt.show()
  
  
  
else:
  plt.xlabel(r'$\nu$')
  plt.ylabel(r'$|I_r(k)|$')
  plt.margins(0.1, 0.1)
  plt.title("Back EMF current spectrum")
  plt.semilogy(freq[:nfft],modfft[:nfft],"b.") # plot modulus of Fourier Transform
  plt.show()

