import matplotlib.pyplot as plt
import numpy as np

def normalizeplot():
	distarray = np.array([0.9006322922498897, 2.7439524066648158, 1.4953113573271992,
		3.14364496729777, 2.580147934431025, 0.5720111882195087, 1.2186095301511557,
		0.6202110860420154, 1.4649405787756613])
	constarray = np.array([.0625, .125, .25, .5, 1.0, 2.0, 4.0, 8.0, 16.0])

	plt.plot(constarray, distarray)
	plt.xlabel("Constant Value")
	plt.ylabel("Average Distance (m) over 1000 iterations")
	plt.title("Avg Distance vs. C-value")
	plt.show()