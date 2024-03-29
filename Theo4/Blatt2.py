from Source import *
# %%
# =============== Aufgabe 1 ===============
from scipy.stats import binom

# binomial dist with p=0.5 and n=100
p = 0.5
n = 100

k = np.arange(n+1)
# k = 50
ways = binom.pmf(k, n, p)

# S/kb
S = np.log(ways)
# E / muB B
E = n - 2 * k

# plot entropy and Energy
plt.plot(k, S, label="S")
plt.plot(k, E, label="E")

# calculate Temperature as derivative of entropy with respect to energy
T = 1/np.gradient(S, E)
# calculate Cv as derivative of energy with respect to temperature
Cv = np.gradient(E, T)
plt.plot(k, T, label="T")
plt.plot(k, Cv, label="Cv")

plt.legend()
plt.xlabel("k")

# %%
# =============== Aufgabe 3 ===============
def S_for_dS(N = 1, V=1):
    S = 3*N/2*np.log(N) + N*np.log(V/N)
    return S


def Stot(ratio, j=0):
    N1 = ratio
    N2 = 1 - ratio
    V1 = ratio
    V2 = 1 - ratio

    if j == 0:
        dS1 = S_for_dS(N1, 1) - S_for_dS(N1, V1)
        dS2 = S_for_dS(N2, 1) - S_for_dS(N2, V2)
        return dS1 + dS2
    else:
        return S_for_dS(N1 + N2, 1) - S_for_dS(N1, V1) - S_for_dS(N2, V2)

y1 = Stot(x)
y2 = Stot(x, 1)
plt.plot(x, y1, label="Ar-He")
plt.plot(x, y2, label="He-He")

plt.xlabel("$x$")
plt.ylabel("$S / k_\mathrm{B}$")
plt.legend()
print(y1/y2)

# %%
# mixing of 3 gases with 2 ratios
def Stot3(ratio1, ratio2):
    N1 = ratio1
    N2 = ratio2
    N3 = 1 - ratio1 - ratio2
    V1 = ratio1
    V2 = ratio2
    V3 = 1 - ratio1 - ratio2

    dS1 = S_for_dS(N1, 1) - S_for_dS(N1, V1)
    dS2 = S_for_dS(N2, 1) - S_for_dS(N2, V2)
    dS3 = S_for_dS(N3, 1) - S_for_dS(N3, V3)
    return dS1 + dS2 + dS3

# ration matrix
x = np.linspace(0.01, 0.99, 100)
y = np.linspace(0.01, 0.99, 100)
X, Y = np.meshgrid(x, y)
Z = Stot3(X, Y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)

ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")
ax.set_zlabel("$S/k_B$")

# project heatmap on xy plane
fig, ax = plt.subplots()
c = ax.pcolormesh(X, Y, Z, cmap="viridis")
fig.colorbar(c, ax=ax)
ax.scatter(0.33, 0.33, color="red", s=100)


plt.show()

# %%
# =============== Aufgabe 4 ===============
from scipy.special import gamma
import math
def Z(N, q, gam=0):
    if gam == 0:
        return gamma(q + N) / (gamma(q + 1) * gamma(N))
    else:
        return math.factorial(q + N - 1) / (math.factorial(q) * math.factorial(N - 1))


fig, ax = plt.subplots(1, 2, figsize=(10, 5))
for (N1, N2, qtot, i) in ((3, 3, 6, 0), (250, 150, 100, 1)):
    if i == 0:
        q = np.arange(qtot + 1, step=0.1)
    else:
        q = np.arange(qtot + 1)
    Ztot = [Z(N1, q, gam=i) * Z(N2, qtot - q, gam=i) for q in q]
    S = np.log(Ztot)
    ax[i].plot(q, S, label=f"N1={N1}, N2={N2}, qtot={qtot}")
    ax[i].axvline(q[np.argmax(S)], color="red", linestyle="--", label=f"$q_1={q[np.argmax(S)]}$")
    ax[i].set_xlabel("$q_1$")
    ax[i].set_ylabel("$S\ (\mathrm{k_B})$")

    print(q[np.argmax(S)])
