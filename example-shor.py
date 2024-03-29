

from quantum import *
n = 15
O1 = Hadamard(4) * Identity(4)
a = 7
L = 2 ** 4

preO2 = [[1 if j==pow(a, i, n) else 0 for i in range(L)] for j in range(L)]
O2 = Identity(4) * operator(preO2)

phi0 = QDit(('|0>', 8))
phi1 = O1(phi0)
phi2 = O2(phi1)
phi3 = QFT(8)(phi2)

obs = observer()
obs.plot(phi3)
