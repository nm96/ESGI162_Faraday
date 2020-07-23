import ode_rk4
import numpy as np
import matplotlib.pyplot as plt

class Motor(ode_rk4.ODE_RK4):
  """A class to compute the trajectory of an electron in an 
     electro-magnetic field non-zero only inside a circle.
  """
  def __init__(self, V0=[], dt=0.1, t0=0):
      """Set the electron ODE parameters
      : param V0     : initial value (as an array or a list) 
      : param dt     : integration time step
      : param t0     : initial time
      : param Ex, Ey : electric field 
      : param B      : magnetic field
      """  
      super().__init__(V0,dt,t0) 

      self.Omega = 2*np.pi*50
      self.Gamma = 0.1
      self.I = 1
      self.R = 1
      self.Phi = 0.1
      self.Kb = 0.001
      self.N = 16
      self.Ns = 6 # unused
      self.Ntrel = 5 # unused
      self.tmax = 15
      self.sf = 1000
      self.eps0 = 0
      self.eps_pre = 0
      self.nu_pre = 0
      self.eps_I = 0
      self.nu_I = 0
      self.eps_G = 0
      self.nu_G = 0
      self.wfG= 0

  def set(self, nu, Gamma, Gamma3, I, R, Phi, Kb, N, Ns, trel, tmax, sf, eps0, eps_pre, nu_pre, eps_I, nu_I, eps_G, nu_G, wfG):
      self.nu = nu
      self.Omega = 2*np.pi*nu
      self.Gamma = Gamma
      self.Gamma3 = Gamma3
      self.I = I
      self.R = R
      self.Phi = Phi
      self.Kb = Kb
      self.N = N
      self.Ns = Ns
      self.trel = trel
      self.tmax = tmax
      self.sf = sf
      self.eps0 = eps0
      self.eps_pre = eps_pre
      self.nu_pre = nu_pre
      self.eps_I = eps_I
      self.nu_I = nu_I
      self.eps_G = eps_G
      self.nu_G = nu_G
      self.wfG = wfG

  def F(self, t, v):
      """ equation to solve: 
         dth/dt = dth
         d dth/dt = sum_i K_i (dth-Omega) sin(theta-Omega t +i pi/N)^2
                   -Gamma dth
      : param t : current time 
      : param v : current function as a vector
      """
      th= v[0]; dth=v[1];

      Phi = self.Phi
      Kt = Phi*Phi/self.R
      omega_nu = 2*np.pi*self.nu_pre
      
      if(self.nu_G < 1e-12): # precession with motor
        angleG = th
      else:
        angleG = self.nu_G*2*np.pi*t

      wave = np.sin(angleG)
      if(self.wfG==1):
        wave = wave**3
      elif(self.wfG==2):
        if(wave>0): wave = 1
        else : wave = -1
      precG = 1+self.eps_G*wave

      dPhi = omega_nu*self.eps_pre*self.Phi*np.sin(omega_nu*t)
      F =-(dth-self.Omega)*Kt*(1+self.eps0+self.eps_pre*np.cos(omega_nu*t))**2\
           *np.sin(th-self.Omega*t)**2\
           + Phi*dPhi/(2*self.R)*np.sin(2*(th-self.Omega*t))       
      for i in range(1,self.N):
        phi_i = i*np.pi/self.N
        dPhi = omega_nu*self.eps_pre*self.Phi*np.sin(omega_nu*t+phi_i)
        F += -(dth-self.Omega)* Kt*\
          (1+self.eps_pre*np.cos(omega_nu*t+phi_i))**2\
          *np.sin(th-self.Omega*t+phi_i)**2\
          + Phi*dPhi/(2*self.R)*np.sin(2*(th-self.Omega*t+phi_i)) 
      F -= self.Gamma*precG*dth+ self.Gamma3*dth*dth*dth

      precI = 1+self.eps_I*np.sin(th-self.nu_I*2*np.pi*t)
      
      eq1 = dth        # equation for d th/dt
      eq2 = F/precI    # equation for d dth/dt
      return(np.array([eq1,eq2]))


  def save(self, fname, tmin, tmax=-1):
      fp = open(fname,"w")
      for i in range(len(self.t_list)):
        if (self.t_list[i] >= tmin) and ((tmax < 0) or (self.t_list[i] <= tmax)):
          t = self.t_list[i]
          th = self.V_list[i][0]
          dth = self.V_list[i][1]
          omega_nu = 2*np.pi*self.nu_pre
          Ir = -self.Kb*np.cos(th)*(\
              self.Phi*(1+self.eps0+self.eps_pre*np.cos(omega_nu*t))\
                  *(dth-self.Omega)*np.sin(th-self.Omega*t)\
                -self.Phi*self.eps_pre*omega_nu*np.sin(omega_nu*t))
          for j in range(1,self.N):
            phi_i = np.pi*j/self.N
            Ir +=  -self.Kb*np.cos(th+phi_i)*(\
              self.Phi*(1+self.eps_pre*np.cos(omega_nu*t+phi_i))\
                  *(dth-self.Omega)*np.sin(th-self.Omega*t+phi_i)\
                -self.Phi*self.eps_pre*omega_nu*np.sin(omega_nu*t+phi_i))

          fp.write("{} {} {} {}\n".format(t, th, dth, Ir))
