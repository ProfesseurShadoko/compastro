


import numpy as np
import pandas as pd


class Galaxy:
    G = 1.0
    
    def __init__(
        self,
        N_particles = 5000,
        galaxy_radius = 1.0,
        galaxy_mass = 1.0,
        galaxy_thickness_over_radius = 0.3 / 15,
        softening_over_radius = 0.01,
        smbh_mass_over_galaxy_mass = 5,
        smbh_radius_over_galaxy_radius = 0.5,
    ):
        self.N = N_particles
        self.R = galaxy_radius
        self.M = galaxy_mass
        self.m = self.M / self.N
        self.z0 = galaxy_thickness_over_radius * self.R
        self.eps = softening_over_radius * self.R
        self.M_smbh = smbh_mass_over_galaxy_mass * self.M
        self.R_smbh = smbh_radius_over_galaxy_radius * self.R
        
        self.rho0 = self.M / (4 * np.pi * self.z0 * self.R**2)
    
    
    def rho(self, x, y, z):
        """
        Returns:
            Exponential disk density profile at a given point (x, y, z)
        """
        r = np.sqrt(x**2 + y**2)
        return self.rho0 * np.exp(-r / self.R) * np.exp(-abs(z) / self.z0)
    
    
    def M_cumulative(self, r, m_smbh = None, r_smbh = None):
        """
        Returns:
            Cumulative mass at a given radius r
        
        If m_smbh and r_smbh are not None, they are taken into account as well
        """
        if m_smbh is None or r_smbh is None:
            return self.M * (1 - (1 + r / self.R) * np.exp(-r / self.R)) # the mass profile without the SMBH
        
        return m_smbh + self.M_cumulative(r) - self.M_cumulative(r_smbh)
    
    
    def v_circular(self, r:np.ndarray, m_cumul:np.ndarray=None):
        """
        Returns:
            circular velocity at a given radius r in order to have a stable orbit
        """
        if m_cumul is None:
            m_cumul = self.M_cumulative(r, self.M_smbh, self.R_smbh) # else, the user might want to provide the cumulative mass from the sample, that is not skewed!
        return np.sqrt(
            (self.G * m_cumul * r**2) / (r**2 + self.eps**2)**(3/2)
        )
        
    
    def p(self, r, z):
        """
        Returns:   
            probability density function of the exponential disk
        """
        norm = 4 * np.pi * self.z0 * self.R**2
        return 1 / norm * self.rho(r, 0, z) / self.rho0
    
    
    def sample_particles(self, n_part:int) -> np.ndarray:
        """
        Returns:
            n_part samples from the exponential disk. No smbh.
        """
        sample = []
        n_sigmas = 10
        max_r = n_sigmas*self.R
        max_z = n_sigmas*self.z0
        while len(sample) < n_part:
            x = np.random.uniform(-max_r, max_r)
            y = np.random.uniform(-max_r, max_r)
            z = np.random.uniform(-max_z, max_z)
            r = np.sqrt(x**2 + y**2)
            
            # let's garantee elipsoidal symmetry
            if (x/max_r)**2 + (y/max_r)**2 + (z/max_z)**2 > 1:
                continue
            
            # let's put no particle in the SMBH
            if r < self.R_smbh:
                continue
            
            if np.random.uniform(0, 1) < self.p(r, z):
                sample.append((x, y, z))
        return np.array(sample)
    
    
    def generate(self) -> pd.DataFrame:
        """
        instatiate self.df with the galaxy particles
        """
        # ----------------- #
        # !-- Positions --! #
        # ----------------- #
        print("Sampling particles...", end=" ", flush=True)
        sample = self.sample_particles(self.N - 1)
        r = np.sqrt(sample[:, 0]**2 + sample[:, 1]**2)
        print("OK.")
        
        # ---------------- #
        # !-- Velocity --! #
        # ---------------- #
        
        print("Computing velocities...", end=" ", flush=True)
        
        # we need to compute experimental m_cumul
        # anyway this is basically m_smbh if m_smbh is high enough
        m_cumul = np.zeros_like(r)
        for i in range(len(r)):
            m_cumul[i] = self.m * np.sum(r < r[i])
        m_cumul += self.M_smbh
        
        sample_velocities = self.v_circular(
            r, # maybe here r3d would be better, but it doesn't change anything in the end, the final issue remains even if z is set systematically to 0.
            m_cumul = m_cumul
        )
        print("OK.")
        
        print("Creating dataframe...", end=" ", flush=True)
        
        df = pd.DataFrame(sample, columns=['x', 'y', 'z'])
        df['v'] = sample_velocities 
        df['r'] = r
        df['m'] = self.m
        df['eps'] = self.eps
        df['phi'] = 0.0
        
        theta = np.arctan2(df['y'], df['x'])
        df['vx'] = -df['v'] * np.sin(theta)
        df['vy'] = df['v'] * np.cos(theta)
        df['vz'] = 0.0
        df['v_theta'] = df['v'] # let's keep v_theta for later

        df['v_circle'] = self.v_circular(df['r'])
        df['v'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
        
        # ------------ #
        # !-- SMBH --! #
        # ------------ #
        
        smbh = [0, 0, 0, 0, 0, self.M_smbh, self.eps, 0, 0, 0, 0, 0, 0] # ['x', 'y', 'z', 'v', 'r', 'm', 'eps', 'phi', 'vx', 'vy', 'vz', 'v_theta', 'v_circle']
        df_smbh = pd.concat([pd.DataFrame([smbh], columns=df.columns), df], ignore_index=True)
        
        df_export = df_smbh[
            ["m", "x", "y", "z", "vx", "vy", "vz", "eps", "phi"]
        ].copy(deep=True)
        
        print("OK.")
        
        self.df = df_export
        return df_export
    
    
    def rotate_galaxy(self, angle:float) -> pd.DataFrame:
        """
        Rotates the galaxy by an angle around 0 (its center) in the x-z plane
        """
        self.df['x'] = self.df['x']*np.cos(angle) - self.df['z']*np.sin(angle)
        self.df['z'] = self.df['x']*np.sin(angle) + self.df['z']*np.cos(angle)

    
    
    def add_system_velocity(self, v_sys:np.ndarray):
        """
        Adds system velocity to the galaxy
        """
        self.df['vx'] += v_sys[0]
        self.df['vy'] += v_sys[1]
        self.df['vz'] += v_sys[2]
    
    def move_galaxy(self, vec:np.ndarray):
        """
        Moves galaxy by a given vector
        """
        self.df['x'] += vec[0]
        self.df['y'] += vec[1]
        self.df['z'] += vec[2]
        
    def get(self) -> pd.DataFrame:
        """
        Returns:
            the galaxy dataframe
        """
        return self.df
    
    def add(self, galaxy:'Galaxy'):
        """
        Adds the df of a second galaxy to the current df
        """
        self.df = pd.concat([self.df, galaxy.get()], ignore_index=True)
    
    def shuffle(self):
        """
        Shuffles df to make it random (to use if the function 'add' was called).
        """
        # shuffle the order of the dataframe
        self.df = self.df.sample(frac=1, replace=False).reset_index(drop=True)
    
    def export(self, filename:str = "files/tests/galaxy/milkyway.txt"):
        """
        Exports the galaxy dataframe to a csv file
        """
        self.df.to_csv(filename, index=True, header=False, sep=" ")
        
    def com(self):
        """
        Shifts the galaxy so the the center of mass is at the origin.
        Modifies the systemic velocity of the galaxy so that the center of mass is at rest.
        """
        # com = pondareated mean of the positions by mass
        com = np.average(self.df[['x', 'y', 'z']], axis=0, weights=self.df['m'])
        self.move_galaxy(-com)
        
        # com_velocity = pondareated mean of the velocities by mass
        com_velocity = np.average(self.df[['vx', 'vy', 'vz']], axis=0, weights=self.df['m'])
        self.add_system_velocity(-com_velocity)
        
        
      

if __name__ == "__main__":
    import os
    while not os.path.exists('files'):
        os.chdir('..')
    print("Current working directory:", os.getcwd())
    
    N = 6000
    smbh_mass = 10 # units of total mass
    smbh_radius = 0.3
    
    galaxy1 = Galaxy(N_particles=N, smbh_mass_over_galaxy_mass = smbh_mass, smbh_radius_over_galaxy_radius=smbh_radius)
    galaxy1.generate()
    
    galaxy2 = Galaxy( # a galaxy half as big
        N_particles = N // 3,
        galaxy_radius = 0.3,
        galaxy_mass = 0.4,
        smbh_mass_over_galaxy_mass=smbh_mass,
        smbh_radius_over_galaxy_radius=smbh_radius
    )
    galaxy2.generate()
    
    # let's put them side by side in the same plane
    r = np.array([5.5, 3.5, 0])
    galaxy2.move_galaxy(r)
    
    v_theta = galaxy1.v_circular(np.linalg.norm(r))
    
    galaxy2.add_system_velocity([-v_theta*2, 0, 0])
    
    galaxy1.add(galaxy2)
    galaxy1.com()
    galaxy1.shuffle()
    galaxy1.export("files/tests/galaxy/milkyway.txt")
    
    
    # os.system("make")
    
        
        