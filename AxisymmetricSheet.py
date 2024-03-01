from math import *
import numpy as np

import matplotlib.animation as animation
class AxisymmetricSheet():

	def __init__(self, Q = 1.2, omega = 0.125, kappa = sqrt(2), xi = 0.5, time_begin = 0, time_end = 1):
		self.omega, self.kappa, self.sigma_0 = omega, omega * kappa, xi * omega /(2*pi)  
		self.k_crit = (self.kappa**2)/(2 * pi * self.sigma_0) 

		self.Q, self.oort_A = Q, self.omega * (1 - 0.25 * (self.kappa/self.omega)**2) 
		
		self.time_step = (2*pi*0.01) * (1/self.kappa)
		self.time = np.arange((2*pi/self.kappa) * time_begin, (2*pi/self.kappa) * time_end, self.time_step)

	## General Functions ##
	## ----------------- ## 

	def __getitem__(self, i):
		return self.density_test_ptle[i, :], self.density_consistent[i,:]

	def save_2_file(self, filename, radii, test = True):
		if test:
			out_array = np.vstack(((1/self.omega + radii), self.density_test_ptle))
			np.savetxt(filename, np.real(out_array), delimiter = ',')
		else:
			out_array = np.vstack(((1/self.omega + radii), self.density_consistent))
			np.savetxt(filename, np.real(out_array), delimiter = ',')


	def inverse_ft(self, array, ks, xs): # This is only in one D (kx), but includes 1/2pi from ky integral 
		transformed, delta_k = np.zeros((array.shape[0], xs.shape[0]), dtype = complex), ks[1]-ks[0] 

		for i, row in enumerate(array):
			for j, x in enumerate(xs):
				transformed[i,j] = (1/(2*pi)) * np.sum(row * np.exp(1j * ks * x)) * delta_k 

		return transformed

	def sigma(self):
		return 0.238# 3.36 * self.Q * self.sigma_0 /self.kappa

	def chi(self, k_x):
		return (self.sigma() * k_x /self.kappa)**2

	def print_time_coeff(self):
		print(f"Time Step: {self.time_step:.5f}\nTime End: {self.time[-1]:.5f}\nNumber of Steps: {self.time.shape[0]}\nSigmaR: {self.sigma()}")

	## Axisymmetric Evolution ## 
	## ---------------------- ## 

	def kernel_axisymmetric(self, t_minus_tp, k_x):
		return (abs(k_x)/self.k_crit) * sin(self.kappa * t_minus_tp) * exp((cos(self.kappa * t_minus_tp) - 1) * self.chi(k_x)) 

	def delta_evolution_axisymmetric(self, k_x, guassian_factor = 1):
		density_test_ptle, density_consistent = np.zeros_like(self.time), np.zeros_like(self.time)
		kernels = np.fromiter(map(self.kernel_axisymmetric, self.time, (k_x for _ in self.time)), dtype = float)

		for time_index in range(1, self.time.shape[0]-1):
			density_test_ptle[time_index] = self.kappa * guassian_factor * kernels[time_index] 

			integrand =  self.time_step * self.kappa * kernels[:time_index] * density_consistent[:time_index]
			integrand[0] *= 0.5

			density_consistent[time_index] = (density_test_ptle[time_index] + np.sum(integrand)) / (1 - 0.5 * self.kappa * kernels[time_index])

		return density_test_ptle[:-1], density_consistent[:-1]


	def guassian_delta_evolution(self, delta_x, x_array): 
		
		guassian_factor = lambda kx : (2*pi)* exp(-0.5 * (kx * delta_x)**2) ## sqrt(2pi) factors cancel from the guassian and delta in ky 
		
		ks = self.k_crit * np.linspace(-10 / delta_x, 10 / delta_x, 500)
		density_test_ptle, density_consistent =  np.zeros((self.time.shape[0]-1 , ks.shape[0])), np.zeros((self.time.shape[0]-1 , ks.shape[0])) 

		for i, kx in enumerate(ks):
			density_test_ptle[:, i], density_consistent[:, i] = self.delta_evolution_axisymmetric(kx, guassian_factor(kx)) 

		self.density_test_ptle, self.density_consistent = self.inverse_ft(density_test_ptle, ks, x_array), self.inverse_ft(density_consistent, ks, x_array)
		#self.density_test_ptle, self.density_consistent = density_test_ptle, density_consistent
		

	## Plotting Functions ##
	## ------------------ ## 
	def density_animation(self, filename = None):
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs1 = plt.subplots(1,1)
		ims = []

		x_array = np.linspace(-3, 3, 500)
		delta = 0.1
		self.guassian_delta_evolution(delta, x_array)

		for i, t, c in zip(self.time, self.density_test_ptle, self.density_consistent):
			itemT, = axs1.plot((x_array+(1/self.omega)), np.imag(t)/(1/self.omega),  color = 'royalblue', animated = True)
			itemC, = axs1.plot((x_array+(1/self.omega)), np.real(t)/((1/self.omega)),  color = 'firebrick', animated = True)
			title =fig.text(.4,.9,(f"Time: {i:.2f}" ))
			ims.append([itemT, itemC, title])


		ani = animation.ArtistAnimation(fig, ims, interval=30)

		if filename == None:
			plt.show()
		else:
			ani.save(filename, writer = writer)
			print(f"Animation saved to: {filename}")




