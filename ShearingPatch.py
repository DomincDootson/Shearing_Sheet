from ShearingSheet import *

class ShearingPatch():

	def __init__(self, k_list, Q, omega = 1, kappa = sqrt(2), xi = 1):
		self.k_list = k_list
		self.shearing_sheets = [ShearingSheet(ky, Q, omega, kappa, xi) for ky in self.k_list]


	def __iter__(self):
		for sheet in self.shearing_sheets:
			yield sheet
	def __getitem__(self, index):
		return self.k_list[index]

	## Toomre Amplification ##
	## -------------------- ## 


	def toomre_amplification_RMS(self, ti = -1.5): # Assume that k = [-k_y, k_y] 

		*_, pos_response = self.shearing_sheets[1].delta_evolution(ti)
		*_, neg_response = self.shearing_sheets[0].delta_evolution(-ti) # To keep same k_x

		top = sqrt((1/pos_response.shape[0]) * np.sum(np.square(pos_response)))
		bot = sqrt((1/neg_response.shape[0]) * np.sum(np.square(neg_response)))

		return np.max(abs(pos_response))/np.max(abs(neg_response))

	def toomre_amplification_Abs(self, ti = -1.5): # Assume that k = [-k_y, k_y] 

		*_, pos_response = self.shearing_sheets[1].delta_evolution(ti)
		*_, neg_response = self.shearing_sheets[0].delta_evolution(-ti) # To keep same k_x

		val1, val2 = np.max(abs(pos_response)), np.max(abs(neg_response))
		return val1, val2 #min(val1, val2), max(val1, val2)

	def toomre_amplification(self):
		ti = np.linspace(-1.5, 0, 10)
		RMS = [self.toomre_amplification_RMS(t) for t in ti]

		return max(RMS)

	## Wave Packet Response ##
	## -------------------- ##

	def wave_perturbation(delta, x_i = 0,  M = 1): ## Assume that it is an impulse 
		# Calculate 
		return 1 


			

	## Response to Cloud ##

