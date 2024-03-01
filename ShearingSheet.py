## This Class needs to be split up into two, one that does indivudal k and one that combines different k 

from math import *
import matplotlib.pyplot as plt
import numpy as np

from Animator import *

class ShearingSheet():
	
	def __init__(self, ky, Q, omega = 1, kappa = sqrt(2), xi = 1, time_begin = 0, time_end = 5):	# ky is in units of k_crit
		
		self.omega, self.kappa, self.sigma_0 = omega, omega * kappa, xi *(omega/(2*pi))
		self.k_crit = (self.kappa**2)/(2 * pi * self.sigma_0) 

		self.ky, self.Q = ky * self.k_crit, Q #ky, Q#
		self.oort_A = self.omega * (1 - 0.25 * (self.kappa/self.omega)**2)

		self.time_step = (2*pi*0.01) * (1/self.kappa)
		self.time = np.arange((2*pi/self.kappa) * time_begin, (2*pi/self.kappa) * time_end, self.time_step)

	## Saving Functions ##
	## ---------------- ##

	def save_2_file(self, filename, array, radii):
		
		out_array = np.vstack(((1/self.omega + radii), array))
		np.savetxt(filename, np.real(out_array), delimiter = ',')
		print(f"Save to: {filename}")

	## Misc ##
	## ---- ##

	def sigma_R(self):
		sigma = 3.36 * self.sigma_0 * self.Q /(self.kappa)
		print(f"Sigma_R: {sigma}")
		return sigma

	def set_Q_with_sigma(self, sigmaR):
		self.Q = self.kappa * sigmaR /(3.36*self.sigma_0)


	## Kernel Functions ## 
	## ---------------- ##

	def S(self, t):
		return sin(self.kappa * t)

	def C(self, t):
		return cos(self.kappa * t)

	def bx(self, t, tp): # Absorb the k/k_crit into def
		return abs(self.ky/self.k_crit) * (self.oort_A * (tp*self.S(tp) - t * self.S(t)) + (self.omega/self.kappa) * (self.C(tp) - self.C(t)))

	def by(self, t, tp):
		return abs(self.ky/self.k_crit) * (self.oort_A * (tp*self.C(tp) - t * self.C(t)) - (self.omega/self.kappa) * (self.S(tp) - self.S(t)))

	def mag_b(self, t, tp):
		return self.bx(t, tp)**2 + self.by(t, tp)**2

	def cx(self, tp):
		return -self.oort_A * tp * self.C(tp) + (self.omega/self.kappa) * (self.S(tp))

	def cy(self, tp):
		return  self.oort_A * tp * self.S(tp) + (self.omega/self.kappa) * (self.C(tp))

	def b_dot_c(self, t, tp):
		return self.bx(t, tp) * self.cx(tp) + self.by(t, tp) * self.cy(tp)

	def kernel(self, t, tp):
		return 4 * self.b_dot_c(t, tp) * exp(-0.572 * self.Q**2 * self.mag_b(t, tp)) / sqrt(1 + 4 * self.oort_A**2 * tp**2)

	def kernel_grid(self, t_array):
		kernels = self.kappa*np.fromiter(map(self.kernel, (t_array[-1] for _ in t_array), t_array), dtype = float)
		kernels[0] *= 0.5
		kernels[-1] *= 0.5
		return kernels


	## Delta Function Evolution ##
	## ------------------------ ##

	def delta_evolution(self, ti): 
		t_end = 5
		t_array = np.arange(ti * pi /self.kappa, t_end * pi /self.kappa, self.time_step)
		density_test_ptle, density_consistent = np.zeros_like(t_array), np.zeros_like(t_array)
		
		for t_index in range(1, t_array.shape[0]-1): 
			kernels = self.kernel_grid(t_array[:t_index+1])
			
			pert = kernels[0] * (2/self.kappa)  
			density_int = np.sum(kernels[:t_index] * density_consistent[:t_index]) * self.time_step

			density_consistent[t_index] = (pert + density_int) / (1 - kernels[-1])
			density_test_ptle[t_index]  = pert

		return (self.kappa/pi) * t_array[:-2], density_test_ptle[:-2],  density_consistent[:-2]

	def delta_amplification_ti(self, ti): # This to return the RMS 
	 response  = self.delta_evolution(ti)
	 return np.amax(response[2])/np.amax(response[1])

	def delta_amplification(self):
		ti = np.linspace(-1.5, -0., 10)
		list_max = [self.delta_amplification_ti(t) for t in ti]
		return max(list_max)


	## Cloud Evolution ## 
	## --------------- ## 

	def mag_k_sq(self, tp):
		return (self.ky**2) * (1 + (2 * self.oort_A * tp)**2)

	def cloud_pert(self, tp, M = 1, delta = 0):
		return  M * exp(-0.5 * self.mag_k_sq(tp) * (delta **2)) 

	def cloud_evolution(self, ti, M = 1, delta = 0): 
		t_end = 5
		t_array = np.arange(ti * pi /self.kappa, t_end * pi /self.kappa, self.time_step)
		density_test_ptle, density_consistent, pert = np.zeros_like(t_array), np.zeros_like(t_array), np.asarray(list(map(self.cloud_pert, t_array, (M for _ in t_array), (delta for _ in t_array))))
		
		for t_index in range(1, t_array.shape[0]-1): 
			kernels = self.kernel_grid(t_array[:t_index+1])
			
			pert_int    = np.sum(kernels * pert[:t_index+1])    * self.time_step
			density_int = np.sum(kernels[:t_index] * density_consistent[:t_index]) * self.time_step

			density_consistent[t_index] = (pert_int + density_int) / (1 - kernels[-1])
			density_test_ptle[t_index]  = pert_int

		return (self.kappa/pi) * t_array[:-2], density_test_ptle[:-2],  density_consistent[:-2]




	## Wave Evolution ##
	## -------------- ##

	def wave_perturbation(self, k_x_0, x_i, delta,  M):
		return pi * np.exp(-np.square(k_x_0*delta)*0.5 - 1j * k_x_0 * x_i)


	def individual_k_wave_evolution(self, x_i, delta, n_points):
		
		k_x_0 = np.linspace(x_i-10*delta, x_i+10*delta, n_points)

		# Set up grids
		t_i, perturbation = k_x_0/(2 * self.oort_A * self.ky), self.wave_perturbation(k_x_0, x_i, delta, 1)
		k_grid, test_k, self_k = np.zeros((self.time.shape[0], n_points), dtype = float), np.zeros((self.time.shape[0], n_points), dtype = complex), np.zeros((self.time.shape[0], n_points), dtype = complex)

		for i, t in enumerate(self.time):
			if i%20==0: 
				print(f"Fraction completed: {t/self.time[-1]:.2f}")
			kernel = self.kappa * np.fromiter(map(self.kernel, t_i +t, t_i), dtype = complex)
			k_grid[i,:] = k_x_0 + 2 * self.oort_A * t * self.ky 
			test_k[i, :] = kernel * perturbation
			
			integral = np.zeros((n_points), dtype = complex)
			for j, t_prime in enumerate(self.time[:i]): 
				integral += self.kappa * self.time_step * self_k[j,:] * np.fromiter(map(self.kernel, t_i + t, t_i + t_prime), dtype = complex)
				if j ==0:
					integral *= 0.5

			self_k[i,:] = (test_k[i,:] + integral)/(1-0.5*self.time_step*self.kappa*np.fromiter(map(self.kernel, t_i +t, t_i+t), dtype = complex))


		return k_grid, test_k, self_k 

	def inverse_ft(self, array, ks, xs): 
		transformed = np.zeros((array.shape[0], xs.shape[0]), dtype = complex)

		for i, row in enumerate(array):
			for j, x in enumerate(xs):
				delta_k = abs(ks[i, 1] - ks[i, 0])
				
				transformed[i,j] = (1/(2*pi)) * np.sum(row * np.exp(1j * ks[i,:] * x)) * delta_k 

		return transformed/(2*pi) # From the inverse FT in the y direction 


	def impulse_density_evolution(self, x_i, delta, n_points = 100, filename = None):
		k_grid, test_k, self_k = self.individual_k_wave_evolution(x_i, delta, n_points)
		
		x_vector = np.linspace(-8*delta, 8*delta, n_points)
		test, cons = 2*np.real(self.inverse_ft(test_k, k_grid, x_vector)), 2*np.real(self.inverse_ft(self_k, k_grid, x_vector)) # Density must be real
		
		if filename == None:
			return x_vector, test, cons  
		else: 
			self.save_2_file(filename+"_Test.csv", test, x_vector)
			self.save_2_file(filename+".csv",      cons, x_vector)
			return x_vector, test, cons  


# animator = Animator(x, [d,c])
# animator.animate()
# sheet = ShearingSheet(-1, 1.2)
# response = sheet.delta_evolution(-1.5)
# plt.plot(response[0], response[2])


