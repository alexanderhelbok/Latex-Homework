from Source import *
from scipy.stats import gaussian_kde


def get_initial_positions(N_atoms, box_size):
    return np.random.uniform(0, box_size, size=(N_atoms, 3))


def get_dist(pos1, pos2):
    return ((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2) ** 0.5


def get_energy(pos):
    en = 0.
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            dist = get_dist(pos[i], pos[j]) ** 4
            en += 4. / dist ** 2 - 4. / dist
    return en


def try_step(pos, iatom, length, box_size):
    radius = np.random.uniform(0, length)
    angle = np.random.uniform(0, 2 * np.pi)
    theta = np.random.uniform(0, np.pi)
    pos[iatom][0] += radius * np.cos(angle) * np.sin(theta)
    pos[iatom][1] += radius * np.sin(angle) * np.sin(theta)
    pos[iatom][2] += radius * np.cos(theta)
    if (pos[iatom][0] > box_size) or (pos[iatom][0] < 0.) or (pos[iatom][1] > box_size) or (pos[iatom][1] < 0.) or (
            pos[iatom][2] > box_size) or (pos[iatom][2] < 0.):
        en = 1e18
    else:
        en = get_energy(pos)
    return (pos, en)


def p(T, N, V):
    return N * T / V


def dreid_simulation(N_atoms, box_size_old, n_steps, temperature):
    global box_size_new, pos_new, en_new
    pos_old = get_initial_positions(N_atoms, box_size_old)
    en_old = get_energy(pos_old)
    for istep in range(n_steps):
        if istep % 1000 == 0:
            print("Step", istep, "Energy:", en_old)
        V_oder_ort = np.random.uniform(0, 1)
        if V_oder_ort < 0.95:  # zu 95 % wird ort eines teilchens geändert, wie bisher
            (pos_new, en_new) = try_step(np.copy(pos_old), np.random.choice(N_atoms), 5., box_size_old)
        if V_oder_ort >= 0.95:  # zu 5 % wird volumen geändert
            box_size_new = box_size_old * np.random.uniform(0.95,1.05)  # volumsänderung um bis zu 5 % mehr oder weniger
            u = np.random.uniform(0, 1)
            for i in range(0, 1, 2):
                pos_new[i] = pos_old[i] * (box_size_new / box_size_old) ** 3
            if (u<np.exp(-((en_new-en_old+p(N_atoms, temperature, box_size_new**3)*(np.abs((box_size_new**3)-box_size_old**3))) / temperature))):
                pos_old = np.copy(pos_new)
                box_size_old = np.copy(box_size_new)
                en_old = en_new
                box_size_old = box_size_new  # übernehmen der neuen größe als alte für nächsten durchlauf
    return box_size_new ** 3, p(N_atoms, temperature, box_size_new ** 3)


temperatures = [1e-5, 2e-5, 3e-5, 0.5, 0.6, 0.78, 1, 200, 300, 350]
pressures = []
volumes = []
for i in temperatures:
    a, b = dreid_simulation(40, 10, 30000, i)
    volumes.append(a)
    pressures.append(b)

# %%
import matplotlib
fig, ax = plt.subplots(figsize=(7, 5))
plt.scatter(pressures, volumes, c=temperatures, cmap="cividis", label="solid", norm=matplotlib.colors.LogNorm(), ec="k", s=50)
# colorbar
# plt.axvline(0.5, color="red", linestyle="--", label="critical point")
# plt.axvline(10 ** -3, color="red", linestyle="--")

plt.colorbar(label="T")
plt.xlabel("$p$")
plt.ylabel("$V$")

plt.xscale('log')
#plt.yscale('log')
plt.tight_layout()
plt.savefig("bilder/phase_diagramm.pdf", bbox_inches="tight")
plt.show()
# def update_plots(ax):
#     ax[0].cla()
#     ax[2].cla()
#
#     ax[0].set_title('Atoms')
#     ax[1].set_title('Energy')
#     ax[2].set_title('Distance distribution')
#
#     ax[1].set_xlabel('Step')
#     ax[1].set_ylabel('Energy')
#
#     ax[2].set_xlabel('Distance')
#     ax[2].set_ylabel('Probability')
#
#     ax[0].set_xlabel('x')
#     ax[0].set_ylabel('y')
#
# # %%
# # MAIN CODE STARTS HERE
#
# N_atoms = 40
# box_size = 10
# n_steps = 30001
# # n_steps = 10001
# # temperature = .1
# # temperature = 10
# temperature = 0.9
# pos_old = get_initial_positions(N_atoms, box_size)
# en_old = get_energy(pos_old)
# snapshot_average = 50
# distances = np.zeros(sum(np.arange(1, N_atoms)) * snapshot_average)
#
# print("Starting energy:", en_old)
# # %%
# fig, ax = plt.subplots(1, 3, figsize=(15, 5))
# plt.ion()
# ax[0].set_xlabel('x')
# ax[0].set_ylabel('y')
# ax[0].set_xlim([-1, 11])
# ax[0].set_ylim([-1, 11])
#
# for istep in range(n_steps):
#     (pos_new, en_new) = try_step(np.copy(pos_old), np.random.choice(N_atoms), 5., box_size)
#     u = np.random.uniform(0, 1)
#     if u < np.exp(-(en_new - en_old) / 1. / temperature):
#         pos_old = np.copy(pos_new)
#         en_old = en_new
#     # calculate distances
#     if istep % snapshot_average == 0:
#         snapshot = 0
#         for i in range(N_atoms):
#             for j in range(i + 1, N_atoms):
#                 distances[snapshot] = get_dist(pos_old[i], pos_old[j])
#                 snapshot += 1
#     # if istep % 10000 == 0:
#     if istep % 100 == 0:
#         update_plots(ax)
#
#         ax[0].scatter(pos_old[:N_atoms, 0], pos_old[:N_atoms, 1], alpha=.2, s=400)
#         ax[1].plot(istep, en_old, 'ro')
#         # create histogram with kde
#         ax[2].hist(distances[distances != 0.], bins=50, alpha=0.5, density=True, color="white", edgecolor='black', histtype='stepfilled', label='Histogram')
#         # plot kde on top
#         kde = gaussian_kde(distances[distances != 0.], bw_method=0.1)
#         x = np.linspace(0, max(distances), 1000)
#         ax[2].plot(x, kde(x), 'r-', label='KDE')
#         ax[0].set_xlim([-1, 11])
#         ax[0].set_ylim([-1, 11])
#         # ax[1].set_yscale = 'log'
#         ax[1].set_ylim(-150, 10)
#         # ax[1].set_ylim(-80, 60)
#         # ax[2].set_ylim(0, .02)
#         plt.tight_layout()
#         plt.pause(0.0001)
#
#         print("Step", istep, "Energy:", en_old)
#
# # %%
# plt.legend()
# plt.savefig("/home/taco/Documents/Latex HÜ/Theo4/bilder/flüssig.pdf", bbox_inches='tight')