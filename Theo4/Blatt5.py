import numpy as np
import matplotlib.pyplot as plt
def get_initial_positions(N_atoms,box_size):
    return np.random.uniform(0,box_size,size=(N_atoms,3)) #changed␣
#   size=(N_atoms,2) to size=(N_atoms,3) for simulation in 3D
def get_dist(pos1,pos2,box_size):
    d = ((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2)**0.5 #␣2D: d = ((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)**0.5
    if (d > box_size): #somehow I sometimes received values > box_size and had to␣insert this fail safe
        return box_size
    else:
        return d
def get_energy(pos,box_size):
    en=0.
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            dist=get_dist(pos[i],pos[j],box_size)**4
            en+=4./dist**2-4./dist
    return en
def try_step(pos,iatom,length,box_size):
    #adapted coordinates for 3D:
    radius = np.random.uniform(0,length)
    angle1 = np.random.uniform(0,2*np.pi) #angle in x-y-plane
    angle2 = np.random.uniform(0,np.pi) #angle measured form z-axis
    pos[iatom][0] += radius*np.sin(angle2)*np.cos(angle1)
    pos[iatom][1] += radius*np.sin(angle2)*np.sin(angle1)
    pos[iatom][2] += radius*np.cos(angle1)
    if ((pos[iatom][0] > box_size) or (pos[iatom][0] < 0.) or (pos[iatom][1] > box_size) or (pos[iatom][1] < 0.) or (pos[iatom][2] > box_size) or (pos[iatom][2] < 0.)):
        en = 1e18
    else:
        en = get_energy(pos, box_size)
    return (pos, en)

    # new function for trying volume change step
def try_step_vol(w, volA, volB, posA, posB):  # A...acceptor, B...donator
    volA_new = volA + volB * w

    volB_new = volB * (1 - w)
    fA = (volA_new / volA) ** (1 / 3)  # factor for scaling coordinates with new volume
    fB = (volB_new / volB) ** (1 / 3)  # factor for scaling coordinates with new volume
    posA_new = posA * fA
    posB_new = posB * fB
    box_sizeA_new = volA_new ** (1 / 3)
    box_sizeB_new = volB_new ** (1 / 3)
    enA_new = get_energy(posA_new, box_sizeA_new)
    enB_new = get_energy(posB_new, box_sizeB_new)
    return (posA_new, posB_new, enA_new, enB_new, box_sizeA_new, box_sizeB_new,volA_new, volB_new)

# new function for trying particle exchange step
def try_step_particle_exchange(posA, posB, N_atomsA, N_atomsB, box_sizeA, box_sizeB):  # system 1 is acceptor
    n = np.random.randint(0, N_atomsB)  # index of particle that changes to system A

    N_atomsA_new = N_atomsA + 1
    N_atomsB_new = N_atomsB - 1
    pos_part = np.random.uniform(0, box_sizeA, size=(1, 3))  # new coordinates of exchange particle in system A
    posA_new = np.append(posA, pos_part, axis=0)  # appending new particle coordinates to positions of system A
    posB_new = np.delete(posB, n, axis=0)  # deleteing particle coordinates from␣ positions of system B
    enA_new = get_energy(posA_new, box_sizeA)
    enB_new = get_energy(posB_new, box_sizeB)
    return (posA_new, posB_new, enA_new, enB_new, N_atomsA_new, N_atomsB_new)

# new function for preparing position and energy parameters for simulation
def prepare_parameters(N_atoms, box_size):
    pos_old = get_initial_positions(N_atoms, box_size)

    en_old = get_energy(pos_old, box_size)
    return (pos_old, en_old)

    # new function for plotting positions of particles
def choose_system(pos1, pos2, en1, en2, N_atoms1, N_atoms2, box_size1, box_size2):
    sys = np.random.choice([0,1]) #random number for determining which system is␣affected by operation or acts as acceptor/donator in the operation
    if sys == 0: #system 1 as acceptor or particle is moved in system 1
        posA = np.copy(pos1)
        posB = np.copy(pos2)
        enA = en1
        enB = en2
        N_atomsA = N_atoms1
        N_atomsB = N_atoms2
        box_sizeA = box_size1
        box_sizeB = box_size2
    else: #system 2 as acceptor
        posA = np.copy(pos2)
        posB = np.copy(pos1)
        enA = en2
        enB = en1
        N_atomsA = N_atoms2
        N_atomsB = N_atoms1
        box_sizeA = box_size2
        box_sizeB = box_size1
    return (posA, posB, enA, enB, N_atomsA, N_atomsB, box_sizeA, box_sizeB, sys)

def MC(n_steps, pos1_old, pos2_old, en1_old, en2_old, N_atoms1, N_atoms2, box_size1, box_size2, temperature):
    N1, N2, V1, V2 = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps)
    print ("Starting energies:", en1_old,en2_old)
    for istep in range(n_steps):
        u = np.random.uniform(0,1) #random number for determining if step gets accepted
        v = np.random.uniform(0,1) #random number for determining which operation will be executed (translation of particle within system, exchange of volume, exchange of particle between systems)
        #with 90% probability particle gets moved in randomly chosen system
        if v <= 0.9:
            posA, _ ,enA , _, N_atomsA, _ , box_sizeA, _, sys = choose_system(pos1_old, pos2_old, en1_old, en2_old, N_atoms1, N_atoms2, box_size1, box_size2)
            (posA_new, enA_new) = try_step(np.copy(posA),np.random.choice(N_atomsA),box_sizeA/2,box_sizeA)
            #checking if step gets accepted
            if(u<np.exp(-(enA_new-enA)/1./temperature)):
            #print("Particle translation step accepted. Acceptor:",sys)
                if sys == 0:
                    pos1_old = posA_new
                    en1_old = enA_new
                else:
                    pos2_old = posA_new
                    en2_old = enA_new
            #with probability of 5% volume gets exchanged with randomly chosen acceptor and donator
        elif (0.90<v<=0.95):
            w = np.random.uniform(1e-3,5e-3) #fraction of volume that is exchanged
            #print('w=',w)
            posA, posB ,enA , enB, N_atomsA, N_atomsB , box_sizeA, box_sizeB, sys = choose_system(pos1_old, pos2_old, en1_old, en2_old, N_atoms1, N_atoms2, box_size1, box_size2)
            volA = box_sizeA**3
            volB = box_sizeB**3
            (posA_new, posB_new, enA_new, enB_new, box_sizeA_new, box_sizeB_new, volA_new, volB_new) = try_step_vol(w,np.copy(volA),np.copy(volB),np.copy(posA),np.copy(posB))
        #checking if step gets accepted
            if (u < np.exp( -(1/temperature) * ((enA_new-enA) + (enB_new-enB)) + N_atomsA * np.log((volA_new-w*volB)/volA_new) + N_atomsB*np.log((volB_new+w*volB)/volB_new))) : #new condition for accepting step
        #print('Volume exchange accepted.')
                if sys == 0:
                    pos1_old = np.copy(posA_new)
                    en1_old = enA_new
                    box_size1 = box_sizeA_new
                    pos2_old = np.copy(posB_new)
                    en2_old = enB_new
                    box_size2 = box_sizeB_new
                else:
                    pos1_old = np.copy(posB_new)
                    en1_old = enB_new
                    box_size1 = box_sizeB_new
                    pos2_old = np.copy(posA_new)
                    en2_old = enA_new
                    box_size2 = box_sizeA_new
        #with probability of 5% randomly chosen particle gets exchanges between systems
        else:
        #print('Particle exchange attempted. Acceptor:',sys)
            posA, posB ,enA , enB, N_atomsA, N_atomsB , box_sizeA, box_sizeB, sys = choose_system(np.copy(pos1_old), np.copy(pos2_old), en1_old, en2_old, N_atoms1, N_atoms2, box_size1, box_size2)
            (posA_new, posB_new, enA_new, enB_new, N_atomsA_new, N_atomsB_new) = try_step_particle_exchange(np.copy(posA), np.copy(posB), N_atomsA, N_atomsB, box_sizeA, box_sizeB)
            volA = box_sizeA**3
            volB = box_sizeB**3
        #checking if step gets accepted
        # if (N_atomsB_new > 0): #don't accept step if donator system has only 1 particle left -> leads to calculations problems
            if (u<np.exp( -(1/temperature) * ((enA_new-enA) + (enB_new-enB)) + np.log(volA * (N_atomsB_new + 1)/(N_atomsA_new * volB) ))) : #new condition for accepting step
            #print('Particle exchange ACCEPTED. Acceptor:',sys)
                if sys == 0:
                    pos1_old = np.copy(posA_new)
                    en1_old = enA_new
                    pos2_old = np.copy(posB_new)
                    en2_old = enB_new
                    box_size2 = box_sizeB
                    N_atoms1 = N_atomsA_new
                    N_atoms2 = N_atomsB_new
                else:
                    pos1_old = np.copy(posB_new)
                    en1_old = enB_new
                    pos2_old = np.copy(posA_new)
                    en2_old = enA_new
                    box_size2 = box_sizeA
                    N_atoms1 = N_atomsB_new
                    N_atoms2 = N_atomsA_new

        if(istep%1000==0):
            print("Step",istep, "N_atoms1:",N_atoms1,"N_atoms2:",N_atoms2,"V1:",box_size1**3,"V2:",box_size2**3,"Energy1:",en1_old,"Energy2:",en2_old)
        #     plot(istep,pos1_old,en1_old,N_atoms1,box_size1,title='System 1',color='blue')
        #     plot(istep,pos2_old,en2_old,N_atoms2,box_size2,title='System 2',color='green')
        V1[istep] = box_size1**3
        V2[istep] = box_size2**3
        N1[istep] = N_atoms1
        N2[istep] = N_atoms2
    return V1, V2, N1, N2



n_steps = 30001
N_atoms1 = 100 #particle number of system 1
N_atoms2 = 10 #particle number of system 2
box_size_ini1 = 10 #box_size of system 1
box_size_ini2 = 100 #box_size of system 2
T = 100 #value for temperature

pos1, en1 = prepare_parameters(N_atoms1, box_size_ini1)
pos2, en2 = prepare_parameters(N_atoms2, box_size_ini2)
V1, V2, N1, N2 = MC(n_steps, pos1, pos2, en1, en2, N_atoms1, N_atoms2, box_size_ini1,box_size_ini2, T) 


# %%
x = np.arange(0, n_steps)
fig, ax = plt.subplots(2, sharex=True)
ax[0].plot(x, V1, c="C0", label="System A")
ax[0].plot(x, V2, c="C1", label="System B")
ax[1].plot(x, N1, c="C0", label="System A")
ax[1].plot(x, N2, c="C1", label="System B")

ax[0].set_ylabel("Volume")
ax[1].set_ylabel("Particles")

plt.legend()
plt.show()
