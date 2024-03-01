from AxisymmetricSheet import * 
from ShearingSheet import * 
from ShearingPatch import * 

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import csv

def save_2_file(filename, data):
	with open(filename, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		
		writer.writerows(data)


def binney_amplification(Q=1.2):
	lambdas = np.linspace(0.1, 5, 50)
	amp = []

	for l in lambdas:
		print(l)
		sheet = ShearingSheet(1/l, Q)
		amp.append(sheet.delta_amplification()) 

	# plt.plot(lambdas, amp)
	# plt.yscale('log', base = 10)
	# plt.show()

	return lambdas, amp

def toomre_amplification(Q = 1.2):
	lambdas = np.linspace(0.1, 5, 50)
	amp = []

	for l in lambdas:
		print(f"Q: {Q}, lambda: {l:.2f}")
		patch = ShearingPatch([-1/l, 1/l], Q, omega = 0.2, xi = 0.5)
		amp.append(patch.toomre_amplification()) 

	return lambdas, amp


def save_amplification_data(filename_B = "Binney_Amplification_Comp.csv", filename_T = "Toomre_Amplification_Comp.csv"):
	Q = [1.3, 1.5]

	toomre_amp = [[*toomre_amplification(q)] for q in Q]
	
	with open(filename_T, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(toomre_amp[0][0])
		for amp in toomre_amp:
			writer.writerow(amp[1])
	print("Now Doing Binney")
	binney_amp = [[*binney_amplification(q)] for q in Q]
	with open(filename_B, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(binney_amp[0][0])
		for amp in binney_amp:
			writer.writerow(amp[1])



def axiysmmetric_response():
	sheet  = AxisymmetricSheet()
	sheet.print_time_coeff()
	x_array = np.linspace(-4, 4, 300)
	delta = 1
	sheet.guassian_delta_evolution(delta, x_array)
	#sheet.save_2_file("../Plotting/sheet_consistant.csv", x_array, test = False)
	print(sheet)

def ky_diff_R():
	l_vec = range(2, 20)
	R_vec = np.linspace(1,10)
	fig, axs = plt.subplots(ncols =2, sharey = True)
	lst1, lst2 = [], []
	for l in l_vec:	
		r =1
		ky = (pi) /(r*sin(pi/l))
		sheet = ShearingSheet(ky, 1.2, omega = 1/r, xi = 0.5)
		lst1.append(ky/sheet.k_crit)

		ky = (l) /(r)
		sheet = ShearingSheet(ky, 1.2, omega = 1/r, xi = 0.5)
		lst2.append(ky/sheet.k_crit)

	axs[0].plot(l_vec, lst1)
	axs[0].plot(l_vec, lst2)
	plt.show()

def max_density_response():
	data = []
	for l in np.linspace(1,8):
		print(l)
		patch = ShearingPatch([l/4, -l/4], 1.5, 0.2, xi = 0.5)
		data.append([l, *patch.toomre_amplification_Abs()])
		
	save_2_file("../Plotting/swing_test_15.csv", data)

#save_amplification_data()


def shearing_sheet_response_ti():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	sheet = ShearingSheet(1, 1)
	t_i = np.linspace(-2, 0.5, 6)
	
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-2., vmax=0.75))
	

	for t in t_i:
		response = sheet.delta_evolution(t)
		label = f"{t:.1f}" if t < 0 else f" {t:.1f}"
		plt.plot(response[0], response[2], label = label, color = cmap.to_rgba(t))
	

	legend = plt.legend(title = r"$\kappa t_{i}/\pi$", title_fontsize =15, fontsize = 15)
	# for t in legend.get_texts():
   	# 	t.set_ha('center')

	plt.xlabel(r"$t\kappa/2\pi$", fontsize = 15)
	plt.ylabel(r"$\tilde{\Sigma}_{1}/\hat{\Sigma}_{e}$", fontsize = 15)
	plt.show()


## Cold disc Amplification ## 
## ----------------------- ##


def cold_disc_amplification(Q, ell):
	v_xi = np.linspace(0.01, 0.50, 40)
	k_over_kcrit = v_xi * ell * 0.5 
	
	amp = []
	for xi, k_norm in zip(v_xi, k_over_kcrit):
		print(f"{Q}, {ell}, {xi}")
		sheet = ShearingSheet(k_norm, Q, xi = xi)
		amp.append(sheet.delta_amplification())
		
	return v_xi, amp

	
def cold_disc_limit(Q = 1.3, filename = "../Plotting/Swing_Data/Cold_Limit/Cold_Limit_Sheet_13.csv"):
	ell = [8, 4, 8]

	with open(filename, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		
		v_xi, amp = cold_disc_amplification(Q, ell[0])
		print(amp)

		writer.writerow(v_xi)
		writer.writerow(amp)

		# for l in ell[1:]:
		# 	v_xi, amp = cold_disc_amplification(Q, l)
		# 	writer.writerow(amp)


max_density_response()
# cold_disc_limit()

#shearing_sheet_response_ti()

# axiysmmetric_response()

# sheet = ShearingSheet(ky = 2.0, Q = 1.288, omega = 0.125, kappa = sqrt(2), xi = 0.5, time_begin = 0, time_end = 2)
# sheet.set_Q_with_sigma(0.238)
# print(sheet.time_step)

# x, d, c = sheet.impulse_density_evolution(0,0.5, filename = "../Plotting/Shearing_Sheet_Comparison_Warm") 







