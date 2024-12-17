import numpy as np
import pandas as pd

N_part = 10_000
N_save = 1000
N_iter = 1000

# !-- Units --! #

G = 1
M_tot = 1 # units of 10^12 M_solar
R = 1
m = M_tot / N_part
z0 = 0.3 / 15 # units of 15kpc
eps = R * 0.01
# !-- Time Scales --! #
t_cross = np.sqrt(R**3 / (G * M_tot))
dt = t_cross / 50
t_relax = N_part / (8 * np.log(R / eps) * t_cross) if eps > 0 else 0


print("Softening Length: ", eps)
print("Crossing Time: ", t_cross)
print("Time step: ", dt)
print("Duration of the simulation: ", N_iter * dt)
print("Relaxation Time: ", t_relax)
print("Thicknes of the disk: ", z0)






rho0 = M_tot / (4 * np.pi * z0 * R**2)

def rho(x, y, z):
    r = np.sqrt(x**2 + y**2)
    return rho0 * np.exp(-r / R) * np.exp(-abs(z) / z0)

def M_cumulative(r):
    return M_tot * (1 - (1 + r / R) * np.exp(-r / R))

def v_circular(r, m_cumul=None):
    if m_cumul is None:
        m_cumul = M_cumulative(r) # else, the user might want to provide the cumulative mass from the sample, that is not skewed!
    return np.sqrt(
        (G * m_cumul * r**2) / (r**2 + eps**2)**(3/2)
    )
    
    
    
def p(r, z = None):
    if z is None:
        return 2 * z0 * p(r, 0)
    
    norm = 4 * np.pi * z0 * R**2
    return 1 / norm * rho(r, 0, z) / rho0


def sample_particles(n_part):
    sample = []
    n_sigmas = 10
    max_r = n_sigmas*R
    max_z = n_sigmas*z0
    while len(sample) < n_part:
        x = np.random.uniform(-max_r, max_r)
        y = np.random.uniform(-max_r, max_r)
        z = np.random.uniform(-max_z, max_z)
        
        # let's garantee elipsoidal symmetry
        if (x/max_r)**2 + (y/max_r)**2 + (z/max_z)**2 > 1:
            continue
        
        r = np.sqrt(x**2 + y**2)
        if np.random.uniform(0, 1) < p(r, z):
            sample.append((x, y, z))
    return np.array(sample)

sample = sample_particles(N_part)




# /!\ CHECK THIS, VELOCITY ISSUES /!\ #

r = np.sqrt(sample[:, 0]**2 + sample[:, 1]**2 + sample[:, 2]**2) # this is 3d r

# we need to redefine m_cumul because the previous one was done on sorted array
m_cumul = np.zeros_like(r)
for i in range(len(r)):
    m_cumul[i] = m * np.sum(r < r[i])

sample_velocities = v_circular( # here 3d r is needed
    r,
    m_cumul = m_cumul
)

df = pd.DataFrame(sample, columns=['x', 'y', 'z'])
df['v'] = sample_velocities
df['r'] = r
df['m'] = m
df['eps'] = eps
df['phi'] = 0.0

# df['v'] is df['v_theta'] # let's compute vx vy vz accrodingly
theta = np.arctan2(df['y'], df['x'])
df['vx'] = -df['v'] * np.sin(theta)
df['vy'] = df['v'] * np.cos(theta)
df['vz'] = 0.0

# let's add to vx vy vz a random component that scales with 0.1 v
#for i in range(len(sample)):
#    df.loc[i, ['vx', 'vy', 'vz']] += np.random.normal(0, 0.1 * df.loc[i, 'v'], 3)

df['v_circle'] = v_circular(df['r'])
df['v'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)

