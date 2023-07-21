# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 22:43:07 2021

@author: Erik
"""

import numpy as np
import matplotlib.pyplot as plt

power_law_data = np.genfromtxt('power_law_data.csv', delimiter = ',')
print(power_law_data)

x = power_law_data[:, 0]
y = power_law_data[:, 1]
error = power_law_data[:, 2]

plt.plot(x, y, 'go')
plt.errorbar(x, y, yerr = error, ecolor = 'red')
plt.title('Power Law Data')
plt.show()

plt.plot(np.log(x), np.log(y), 'go')
plt.errorbar(np.log(x), np.log(y), yerr = np.log(error), ecolor = 'red')
plt.title('Power Law Data Log-Log Form')
plt.show()