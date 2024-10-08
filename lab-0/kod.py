import matplotlib.pyplot as plt
import numpy as np

wartosci = np.linspace(0,10,100)
plt.plot(wartosci, np.sin(wartosci))
plt.savefig("f")
plt.show()