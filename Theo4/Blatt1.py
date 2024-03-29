from Source import *
#%%
# =============== Aufgabe 1 ===============
# 2d random walk with 10000 steps and 10000 runs
n_steps = 10000
n_runs = 10000
step_size = 1

# create a random pair of numbers between -1 and 1 for x and y direction
x = np.random.choice([-1, 1], (n_runs, n_steps))
y = np.random.choice([-1, 1], (n_runs, n_steps))
# create a mask for the random walk
mask = np.random.choice([0, 1], (n_runs, n_steps))

# merge x with mask and y with !mask
x = np.where(mask, 0, x)
y = np.where(mask == 0, 0, y)

# calculate the cumulative sum of the random walk
xcum = np.cumsum(x, axis=1)
ycum = np.cumsum(y, axis=1)
# %%
cm = plt.get_cmap("viridis")
# plot the random walk in 2d
# use color gradient to show the time evolutionu using scatter plot
fig, ax = plt.subplots(figsize=(6, 5))
ax.set_prop_cycle(color=[cm(1.*i/(n_steps-1)) for i in range(n_steps-1)])
for i in range(n_steps-1):
    ax.plot(xcum[0][i:i+2], ycum[0][i:i+2])

# mark start and end point
ax.scatter(xcum[0][0], ycum[0][0], label="start", c="C3", zorder=5)
ax.scatter(xcum[0][-1], ycum[0][-1], label="end", c="C0", zorder=5, marker="x")

# add dashed lines at x=0 and y=0
ax.axvline(0, color="gray", linestyle="--", zorder=0)
ax.axhline(0, color="gray", linestyle="--", zorder=0)

# set square limits for the plot
# set_square_lim(ax)
# set square aspect

# add colorbar to show time evolution
sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(vmin=0, vmax=n_steps-1))
sm.set_array([])
fig.colorbar(sm, ax=ax, label="time step")

plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
# plt.savefig("bilder/random_walk_2d.pdf", bbox_inches="tight")
plt.show()
# %%
# construct 2d histogram of the end points of the random walks and compare it to the analytical solution
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(10, 5))
h = ax1.hist2d(xcum[:, -1], ycum[:, -1], bins=100)
# plot 2d normal distribution with mean 0 and variance n_steps using contourf plot
x = np.linspace(-n_steps, n_steps, 1000)
y = np.linspace(-n_steps, n_steps, 1000)
X, Y = np.meshgrid(x, y)
Z = np.exp(-1/2 * (X**2 + Y**2) / n_steps) / (2 * np.pi * n_steps)
ax2.contourf(X, Y, Z, cmap="viridis", levels=100)
# ax2.hist2d(xcum[:, -1], ycum[:, -1], bins=100, norm=LogNorm())
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
# plt.savefig("bilder/random_walk_2d_hist.pdf", bbox_inches="tight")
plt.show()

# %%
# use seabron to create scatterplot of endpoint with marginal histograms
df = pd.DataFrame({"x": xcum[:, -1], "y": ycum[:, -1]})
g = sns.jointplot(data=df, x="x", y="y", kind="hex", marginal_kws=dict(bins=100, stat='density'))

# set square limits for the joint plot
max_lim = set_square_lim(g.ax_joint)

X = np.linspace(-max_lim, max_lim, 1000)

ax_marg_x = g.ax_marg_x
ax_marg_y = g.ax_marg_y

# plot the marginal histograms and the analytical normal distribution
ax_marg_x.plot(X, norm.pdf(X, 0, np.sqrt(n_runs/2)), color='red')
ax_marg_y.plot(norm.pdf(X, 0, np.sqrt(n_runs/2)), X, color='red')

# remove tick labels on marginal axes (and main)
ax_marg_x.tick_params(bottom=False, which='both', top=False)
ax_marg_y.tick_params(left=False, which='both', right=False)
g.ax_joint.tick_params(top=False, right=False, which='both')


print(np.mean(df.x), np.std(df.x), np.mean(df.y), np.std(df.y))
plt.tight_layout()
# plt.savefig("bilder/random_walk_2d_joint.pdf", bbox_inches="tight")
plt.show()

# %%
from scipy.stats import binom
# =============== Aufgabe 2 ===============
# binomial dist with p=0.5 and n=100
p = 0.5

for n in [100, 1000]:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # calculate amount of ways to get k heads
    k = np.arange(n+1)
    ways = binom.pmf(k, n, p)

    # calculate std dev of binomial distribution assuming distribution is normal
    std_dev = np.sqrt(n * p * (1 - p))
    print(f"n={n}: std_dev = {std_dev}")

    # plot the binomial distribution
    ax1.bar(k, ways, edgecolor="black")
    ax1.text(0.2, 0.8, f"$\\sigma ={std_dev:.2f}$", transform=ax1.transAxes, fontsize=14, ha="center")
    ax1.set_xlabel("k")
    ax1.set_ylabel("P(k)")

    if n == 100:
        ax1.set_xlim(30, 70)
    else:
        ax1.set_xlim(450, 550)

    # plot entropy using S = ln(W) with W being the amount of ways to get k heads
    ax2.plot(k, np.log(ways), label=f"n={n}")
    ax2.set_xlabel("k")
    ax2.set_ylabel("$S / k_B$")

    plt.tight_layout()
    # plt.savefig(f"bilder/binomialn{n}.pdf", bbox_inches="tight")
    plt.show()

# %%
# =============== Aufgabe 3 ===============
from scipy.stats import poisson, norm, binom
# create a 20x20 matrix filled with n epsilon = k
dim = 20
dim3 = 2
k = 2
n_iterations = 5000
n_runs = 1000
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# create array to hold the amount of cells with the same energy
tot_counts = np.zeros((int(np.sqrt(n_iterations)),))

for run in range(n_runs):
    matrix = np.full((dim, dim), k)
    iter = 0
    while iter < n_iterations:
        # choose a random point in the matrix
        i, j = np.random.randint(0, dim, 2)
        # choose a random direction
        direction = directions[np.random.choice(4)]
        # calculate the new energy
        matrix[i, j] -= 1
        matrix[(i+direction[0]) % dim, (j+direction[1]) % dim] += 1
        if matrix[i, j] < 0:
            matrix[i, j] += 1
            matrix[(i+direction[0]) % 20, (j+direction[1]) % 20] -= 1
            iter -= 1
        iter += 1

    # print(matrix)
    # count cells with same energy
    unique, counts = np.unique(matrix, return_counts=True)
    # print(dict(zip(unique, counts)))
    tot_counts[unique] += counts

print(tot_counts)

# drop trailing 0 from tot_counts
firstzero = np.where(tot_counts > 0)[0][-1] + 1
tot_counts = tot_counts[:firstzero]
rel_counts = tot_counts/np.sum(tot_counts)
# %%
# create an array of matrices
arr = np.full((2, 2, 2), k)
# fill the array with random numbers
arr = np.random.choice([1, 2, 3, 4, 5], (3, 2, 2))
print(arr)
# %%
rand = np.random.randint(0, dim, size=(2, 3))
indexs = np.arange(0, 3), rand[0], [0, 0, 0]
# arr[indexs] = 0
# print(arr[indexs])
# %%
from astropy.constants import k_B
from astropy import units as u
X, x2 = np.arange(0, len(tot_counts)), np.arange(0, len(tot_counts), 0.2)

def exponential(x, a, b):
    return a * np.exp(-x / b)

popt, pcov = curve_fit(exponential, X, rel_counts)

T = popt[1]/k_B * 2 * u.eV
# convert T to si units
T = T.to(u.K)
print(f"T = {T:.2f}")

# plot the amounts as pdf
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(X, rel_counts, edgecolor="black", alpha=0.6, zorder=0, width=1)
# ax.plot(unique, norm.pdf(unique, k, np.sqrt(n_iterations)/(dim/np.sqrt(dim3))), label="Gaussverteilung", c="C2")
ax.scatter(X, exponential(X, *popt), c="gold", edgecolor="black", s=50)
ax.plot(x2, exponential(x2, *popt), c="black", zorder=0, ls="--")

ax.text(0.5, 0.5, f"$P(E) = {popt[0]:.2f} \cdot e^{{-{popt[1]:.2f} \cdot E}}$", transform=ax.transAxes)

ax.set_xlabel("E")
ax.set_ylabel("P(E)")
ax.set_xticks(np.arange(0, len(tot_counts), 2,  dtype=int))

plt.tight_layout()
# plt.legend()
# plt.savefig("bilder/energy_distribution.pdf", bbox_inches="tight")
plt.show()



