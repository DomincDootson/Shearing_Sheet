import matplotlib.pyplot as plt 
import matplotlib.animation as animation

class Animator(object):
	def __init__(self, x_values, y_values, time_step = 1):
		
		self.x_values = x_values
		self.y_values = y_values

		self.n_time, self.time_step = y_values[0].shape[0], time_step


	def animate(self, filename = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		for time in range(self.n_time):
			lst = []
			for i, array in enumerate(self.y_values):
				color = "firebrick" if (i==0) else "royalblue"
				line, = axs.plot(self.x_values, array[time,:], animated = True, color = color)
				lst.append(line)
			

			lst.append(fig.text(.4,.9,(f"Time: {time * self.time_step:.2f}")))
			ims.append(lst) 


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			print("Animation saved to: " + filename)
			ani.save(filename, writer = writer)
		else:
			plt.show()
		