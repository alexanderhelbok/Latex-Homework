import uncertainties as unc
import scipy.constants as const

e = const.e
epsilon = const.epsilon_0
pi = const.pi
R0 = unc.ufloat(1.3, 0.1)*1e-15

ac = 3*e**2/(20*pi*epsilon*R0)
print(f"{ac:.1uS} J")   # J
print(f"{ac/e:.1uS} eV") # eV
