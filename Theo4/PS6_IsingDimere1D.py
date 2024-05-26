 #Theo4, FF.Locker, ising model 1D
import numpy as np
import matplotlib.pyplot as plt


#Evaluate Probability Distribution Function exp(-H/T), kb absorbed into T.
def get_probability(delta_energy, Temperature):
    return np.exp(-delta_energy / Temperature)

#simple calculation of the mean of a 1d field quantity
def mean(field):
    return(sum(field)/len(field))

#calculate the total energy by summing up the pair energies
def get_energy(spins):
    energy=0
    for i in range(len(spins)-1):
        if i%2==0:
            energy=energy-interaction*spins[i]*spins[i+1]
    return energy


def delta_energy(spins,random_spin):
    #If you do flip one random spin position, the change in energy is:
    if random_spin%2==0: # i can update both, but to stay in the dimer one has to trat odd and even ones. No periodic BC.
        old = -interaction*(spins[random_spin]*spins[random_spin+1]) 
        new = -interaction*((-1)*spins[random_spin]*spins[random_spin+1])
    else:
        old = -interaction*(spins[random_spin-1]*spins[random_spin]) 
        new = -interaction*((-1)*spins[random_spin-1]*spins[random_spin])
    return new-old

def metropolis(spins,L, MC_samples, Temperature, interaction):
    dummy=np.empty([L])
     
    Beta = Temperature**(1)
    data = []
    magnetization=[]
    energy=[]
    #Introducing Metropolis Hastings Algorithim
    
    for i in range(MC_samples):
        #Each Monte Carlo step consists in L random spin moves
        for j in range(L):
            #Choosing a random spin (watch out for even and odd ones in energy of dimer)
            random_spin=np.random.randint(0,L,size=(1))
            #Compuing the change in energy of this spin flip
            delta=delta_energy(spins,random_spin)

            #Metropolis accept-rejection:
            if delta<0:
                #Accept the move if its negative
                spins[random_spin]=-spins[random_spin]
                #print('change')
            else:
                #If its positive, we compute the probability
                probability=get_probability(delta,Temperature)
                random=np.random.rand()
                if random<=probability:
                    #Accept the move
                    spins[random_spin]=-spins[random_spin]

        data.append(list(spins))

        #Afer the MC step, we measure the system
        magnetization.append(sum(spins))
        energy.append(get_energy(spins))
    return data,magnetization,energy




# setting up problem
L = 100 # size of system, choose even to give evryone a friend :))
MC_samples = 100 # number of samples
Temperature = 1 # "temperature" parameter
interaction = 1 # Strength of interaction between nearest neighbours + 1:ferromagnetische interaktion, -1:antiferro.




# running MCMC
# intializing
#Spin Configuration
    
spins = np.random.choice([-1,1],L)

data = metropolis(spins,L, MC_samples, Temperature, interaction)

# Plotting
plt.figure(figsize=(10,10))

plt.subplot(3,1,1)
cbar = plt.imshow(np.transpose(data[0]),vmax=1, vmin=-1)
plt.colorbar(cbar,ticks=[-1,1])
plt.axis('tight')
plt.ylabel('Space',fontdict={'size':20})
plt.title('Evolution of a 1-D Ising Dimer model',fontdict={'size':20})

m0=mean(data[1])


plt.subplot(3,1,2)
plt.plot(data[1],'r')
plt.xlim((0,MC_samples))
plt.ylabel('Magnetisation, mean=%d' %m0,fontdict={'size':20})
plt.xlabel('MC Update',fontdict={'size':10})

plt.subplot(3,1,3)
plt.plot(data[2],'r')
plt.xlim((0,MC_samples))
plt.ylabel('Energy',fontdict={'size':20})


plt.show()