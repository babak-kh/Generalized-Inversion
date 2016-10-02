import numpy as np
import matplotlib.pyplot as plt

magnitudes = np.genfromtxt('/home/babak/PycharmProjects/Simulation/GeneralInversion/Inputs/2046-L.txt')
plt.plot(magnitudes)
plt.show()