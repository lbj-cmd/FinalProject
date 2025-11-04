#Cassandra Bodin
#Phys305
#Final Project: Asteroid Version

#run by "python FinalProjectCode.py"

#Goal of Project:
   #Recreate the threebody problem we did for homework, for Jupiter the Sun and asteroids.
   #Instead of just 4 asteroids though I want to add a large amount of ateroids, maybe 100-1000, all with different semi major axes lengths
   #Plot and animate each in different colors
   #Add code to take the data for the asteroids and use it to plot a histogram showing the Kirkwood Gaps

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec

G = 6.674e-11
GMs = 4*(math.pi**2)
m_sun = 1.989e30

# Jupiter parameters
ecc_jup = 0.0489
semimajor_jup = 5.2044
a_jup = semimajor_jup
perihelion_jup = a_jup*(1-ecc_jup)
m_jup = 1.8982e27
vp_jup = math.sqrt(GMs) * math.sqrt((1+ecc_jup) / (a_jup*(1-ecc_jup)) * (1+(m_jup/m_sun)))

# Calculate MMR positions for Jupiter
# Formula: a = a_jup * (p/q)^(2/3)
def calculate_mmr_positions():
    mmr_ratios = [(3, 1), (5, 2), (7, 3), (2, 1)]  # p:q ratios
    mmr_positions = []
    for p, q in mmr_ratios:
        a = a_jup * (p/q)**(2/3)
        mmr_positions.append((p, q, a))
    return mmr_positions

# Fast Lyapunov Indicator (FLI) calculation
def calculate_fli(a, e, t_integration=100, h=1e-2):
    """
    Calculate FLI for a test particle with given semi-major axis and eccentricity
    """
    # Initial conditions for the test particle
    perihelion = a * (1 - e)
    vp = math.sqrt(GMs) * math.sqrt((1 + e) / (a * (1 - e)))
    x = -perihelion
    y = 0
    vx = 0
    vy = vp
    
    # Initial conditions for a nearby particle (1e-10 AU difference in x)
    x_near = x + 1e-10
    y_near = y
    vx_near = vx
    vy_near = vy
    
    # Jupiter initial conditions
    vx2 = 0
    vy2 = vp_jup
    x2 = -perihelion_jup
    y2 = 0
    
    t = 0
    max_fli = 0
    
    while t < t_integration:
        # Calculate distances
        r2 = math.sqrt(x2**2 + y2**2)
        r = math.sqrt(x**2 + y**2)
        rrel = math.sqrt((x - x2)**2 + (y - y2)**2)
        
        r_near = math.sqrt(x_near**2 + y_near**2)
        rrel_near = math.sqrt((x_near - x2)**2 + (y_near - y2)**2)
        
        # Update velocities for test particle
        vx = vx + h * (-GMs * x / r**3 - GMs * (m_jup/m_sun) * (x - x2) / rrel**3)
        vy = vy + h * (-GMs * y / r**3 - GMs * (m_jup/m_sun) * (y - y2) / rrel**3)
        
        # Update velocities for nearby particle
        vx_near = vx_near + h * (-GMs * x_near / r_near**3 - GMs * (m_jup/m_sun) * (x_near - x2) / rrel_near**3)
        vy_near = vy_near + h * (-GMs * y_near / r_near**3 - GMs * (m_jup/m_sun) * (y_near - y2) / rrel_near**3)
        
        # Update velocities for Jupiter
        vx2 = vx2 + h * (-GMs * x2 / r2**3)
        vy2 = vy2 + h * (-GMs * y2 / r2**3)
        
        # Update positions for test particle
        x = x + vx * h
        y = y + vy * h
        
        # Update positions for nearby particle
        x_near = x_near + vx_near * h
        y_near = y_near + vy_near * h
        
        # Update positions for Jupiter
        x2 = x2 + vx2 * h
        y2 = y2 + vy2 * h
        
        # Calculate distance between test particle and nearby particle
        delta_x = x_near - x
        delta_y = y_near - y
        delta = math.sqrt(delta_x**2 + delta_y**2)
        
        # Normalize the nearby particle position to keep delta small
        if delta > 1e-5:
            factor = 1e-10 / delta
            x_near = x + delta_x * factor
            y_near = y + delta_y * factor
            vx_near = vx + (vx_near - vx) * factor
            vy_near = vy + (vy_near - vy) * factor
            delta = 1e-10
        
        # Update FLI
        if delta > 0:
            current_fli = math.log10(delta / 1e-10)
            if current_fli > max_fli:
                max_fli = current_fli
        
        t += h
    
    return max_fli

# Precompute chaos map on (a, e) grid
def precompute_chaos_map(a_min=2.0, a_max=3.6, e_min=0.0, e_max=0.4, grid_size=50):
    """
    Precompute FLI values for a grid of (a, e) points
    """
    a_grid = np.linspace(a_min, a_max, grid_size)
    e_grid = np.linspace(e_min, e_max, grid_size)
    fli_grid = np.zeros((grid_size, grid_size))
    
    print("Precomputing chaos map...")
    for i in range(grid_size):
        for j in range(grid_size):
            a = a_grid[i]
            e = e_grid[j]
            fli = calculate_fli(a, e)
            fli_grid[i, j] = fli
        print(f"Progress: {i+1}/{grid_size}")
    
    return a_grid, e_grid, fli_grid

# Generate real asteroids with initial positions and velocities
def generate_asteroids(n_asteroids=100):
    """
    Generate real asteroids with random semi-major axes and eccentricities,
    and calculate initial positions and velocities
    """
    semimajor_a = 1.5 * np.random.random(n_asteroids) + 2.0  # between 2 and 3.5 AU
    ecc_a = 0.3 * np.random.random(n_asteroids)  # between 0 and 0.3
    
    x = np.array([])
    y = np.array([])
    vx = np.array([])
    vy = np.array([])
    
    m_a = 1e10  # mass of asteroids (arbitrary since << Msun)
    
    for i in range(n_asteroids):
        a_a = semimajor_a[i]
        e_a = ecc_a[i]
        perihelion = a_a * (1 - e_a)
        
        # Initial parameters for asteroid (perihelion at x=-perihelion, y=0)
        vp = math.sqrt(GMs) * math.sqrt((1 + e_a) / (a_a * (1 - e_a)) * (1 + (m_a / m_sun)))
        astr_vx = 0
        astr_vy = vp
        astr_x = -perihelion
        astr_y = 0
        
        x = np.append(x, astr_x)
        y = np.append(y, astr_y)
        vx = np.append(vx, astr_vx)
        vy = np.append(vy, astr_vy)
    
    return semimajor_a, ecc_a, x, y, vx, vy

def main():
    """
    Main function to run the real-time resonance and chaos dynamics map
    """
    # Precompute chaos map
    a_grid, e_grid, fli_grid = precompute_chaos_map()
    
    # Generate real asteroids with initial conditions
    semimajor_a, ecc_a, x, y, vx, vy = generate_asteroids(n_asteroids=100)
    
    # Calculate MMR positions
    mmr_positions = calculate_mmr_positions()
    
    # Simulation parameters
    h = 1e-2
    t = 0
    t_final = 200
    count = 0
    
    # Jupiter's initial conditions
    vx2 = 0
    vy2 = vp_jup
    x2 = -perihelion_jup
    y2 = 0
    
    # Create UI with two side-by-side windows
    fig = plt.figure(figsize=(15, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    
    # Left: Real-time histogram
    ax_hist = plt.subplot(gs[0])
    ax_hist.set_title('Real-time Asteroid Semi-major Axis Distribution')
    ax_hist.set_xlabel('Semi-major Axis (AU)')
    ax_hist.set_ylabel('Number of Asteroids')
    ax_hist.set_xlim(2.0, 3.6)
    
    # Plot MMR positions as vertical dashed lines
    for p, q, a in mmr_positions:
        ax_hist.axvline(x=a, color='r', linestyle='--', alpha=0.7)
        ax_hist.text(a + 0.01, 0.95, f'{p}:{q}', transform=ax_hist.get_xaxis_transform(),
                    rotation=90, va='top', ha='left', color='r')
    
    # Right: Chaos map
    ax_chaos = plt.subplot(gs[1])
    ax_chaos.set_title('Chaos Map (a vs e)')
    ax_chaos.set_xlabel('Semi-major Axis (AU)')
    ax_chaos.set_ylabel('Eccentricity')
    ax_chaos.set_xlim(2.0, 3.6)
    ax_chaos.set_ylim(0.0, 0.4)
    
    # Plot chaos map
    im = ax_chaos.pcolormesh(a_grid, e_grid, fli_grid.T, cmap='hot', shading='auto')
    cbar = plt.colorbar(im, ax=ax_chaos)
    cbar.set_label('Fast Lyapunov Indicator (FLI)')
    
    # Overlay real asteroids on chaos map
    scat = ax_chaos.scatter(semimajor_a, ecc_a, color='blue', s=10, alpha=0.8)
    
    # Update histogram and chaos map in real-time
    def update_hist(frame):
        nonlocal t, count, x, y, vx, vy, x2, y2
        
        current_semimajor_a = np.copy(semimajor_a)
        current_ecc_a = np.copy(ecc_a)
        
        if t < t_final:
            # Run one time step of the simulation
            
            # Calculate distances for Jupiter
            r2 = math.sqrt(x2**2 + y2**2)
            
            # Update velocities for asteroids
            for i in range(len(x)):
                r = math.sqrt(x[i]**2 + y[i]**2)
                rrel = math.sqrt((x[i] - x2)**2 + (y[i] - y2)**2)
                
                vx[i] += h * (-GMs * x[i] / r**3 - GMs * (m_jup/m_sun) * (x[i] - x2) / rrel**3)
                vy[i] += h * (-GMs * y[i] / r**3 - GMs * (m_jup/m_sun) * (y[i] - y2) / rrel**3)
            
            # Update velocities for Jupiter
            vx2 += h * (-GMs * x2 / r2**3)
            vy2 += h * (-GMs * y2 / r2**3)
            
            # Update positions for asteroids
            for i in range(len(x)):
                x[i] += vx[i] * h
                y[i] += vy[i] * h
            
            # Update positions for Jupiter
            x2 += vx2 * h
            y2 += vy2 * h
            
            # Update time and count
            t += h
            count += 1
            
            # Calculate current semimajor axis and eccentricity for each asteroid
            for i in range(len(x)):
                # Calculate orbital energy E
                v = math.sqrt(vx[i]**2 + vy[i]**2)
                r = math.sqrt(x[i]**2 + y[i]**2)
                E = 0.5 * v**2 - GMs / r
                
                # Calculate semimajor axis a
                a = -GMs / (2 * E)
                current_semimajor_a[i] = a
                
                # Calculate angular momentum h
                h_orb = x[i] * vy[i] - y[i] * vx[i]
                
                # Calculate eccentricity e
                e = math.sqrt(1 + (2 * E * h_orb**2) / (GMs**2))
                current_ecc_a[i] = e
        
        # Update histogram
        ax_hist.clear()
        ax_hist.set_title('Real-time Asteroid Semi-major Axis Distribution')
        ax_hist.set_xlabel('Semi-major Axis (AU)')
        ax_hist.set_ylabel('Number of Asteroids')
        ax_hist.set_xlim(2.0, 3.6)
        
        # Plot histogram with current a values
        ax_hist.hist(current_semimajor_a, bins=50, histtype='bar', color='b', alpha=0.7)
        
        # Replot MMR positions
        for p, q, a in mmr_positions:
            ax_hist.axvline(x=a, color='r', linestyle='--', alpha=0.7)
            ax_hist.text(a + 0.01, 0.95, f'{p}:{q}', transform=ax_hist.get_xaxis_transform(),
                        rotation=90, va='top', ha='left', color='r')
        
        # Update chaos map title with current time
        ax_chaos.set_title(f'Chaos Map (a vs e) - Time: {t:.2f} years')
        
        # Update scatter plot on chaos map
        scat.set_offsets(np.c_[current_semimajor_a, current_ecc_a])
        
        return ax_hist, scat
    
    # Animation
    ani = FuncAnimation(fig, update_hist, frames=100, interval=1000, blit=False)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

def Asteroids():
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation


    G = 6.674e-11
    GMs = 4*(math.pi**2)
    m_sun = 1.989e30

#info Jupiter
    ecc_jup = 0.0489
    semimajor_jup = 5.2044
    a_jup = semimajor_jup
    perihelion_jup = a_jup*(1-ecc_jup)
    m_jup = 1.8982e27
    vp_jup = math.sqrt(GMs) * math.sqrt((1+ecc_jup) / (a_jup*(1-ecc_jup)) * (1+(m_jup/m_sun)))


#info Asteroids
    #generating the asteroids
    n=np.array([])
    for i in range(0,10): #change end value depending on number of atseroids you want to create
        number= 1.5*np.random.random()+2.0 #generate random numbers between 2 and 3.5
        i+=1
        n=np.append(n,number)

    m_a1 = 1e10 #mass of the asteroids are arbitarary since <<Msun

    x=np.array([])
    y=np.array([])
    vx=np.array([])
    vy=np.array([])
    semimajor_a=np.array([])

    for i in range(n.size):
        semimajor_a1 = n[i]
        a_a1 = semimajor_a1
        ecc_a1 = 0 #circular
        perihelion_a1 = a_a1*(1-ecc_a1)


       #initial parameters for asteroid
        vp_a1 = math.sqrt(GMs) * math.sqrt((1+ecc_a1) / (a_a1*(1-ecc_a1)) * (1+(m_a1/m_sun)))
        astr_vx = 0
        astr_vy = vp_a1
        astr_x = -perihelion_a1
        astr_y = 0

        x=np.append(x,astr_x)
        y=np.append(y,astr_y)
        vx=np.append(vx,astr_vx)
        vy=np.append(vy,astr_vy)
        semimajor_a=np.append(semimajor_a,semimajor_a1)

   #initial parameters for Jupiter
    vx2 = 0
    vy2 = vp_jup
    x2 = -perihelion_jup
    y2 = 0

   # steps and time (in years)
    h = 1e-2
    t_final = 200 #need a ton of orbits to see what is happening
    t_initial = 0
    t = t_initial
    count = 0

   #variables for graphing asteroid and variable dependent variables
    xarr = np.array([])
    yarr = np.array([])
    xarr = np.append(xarr,x)
    yarr = np.append(yarr,y)
    # r = math.sqrt(x**2 + y**2)
    # v = math.sqrt(vx**2 + vy**2)

   #variables for graphing jupiter and variable dependent variables
    xarr2 = np.array([])
    yarr2 = np.array([])
    xarr2 = np.append(xarr2,x2)
    yarr2 = np.append(yarr2,y2)
    r2 = math.sqrt(x2**2 + y2**2)
    v2 = math.sqrt(vx2**2 + vy2**2)

    while (t < t_final):
#radial values
        r2 = math.sqrt(x2**2 + y2**2)
        for i in range(n.size):
            r = math.sqrt(x[i]**2 + y[i]**2)
            rrel = math.sqrt( (x[i]-x2)**2 + (y[i]-y2)**2)

#velocity values
# vx, vy is asteroid
# vx2, vy2 is jupiter
            vx[i] = vx[i] +h*(-GMs*x[i]/r**3- GMs*(m_jup/m_sun)*(x[i]-x2)/rrel**3)
            vy[i] = vy[i]+h*(-GMs*y[i]/r**3- GMs*(m_jup/m_sun)*(y[i]-y2)/rrel**3)
        vx2 = vx2 +h*(-GMs*x2/r2**3)#- GMs*(m_jup/m_sun)*(x2-x)/r2**3) #add jupiter mass dependence on asteroid by uncommenting
        vy2 = vy2 +h*(-GMs*y2/r2**3)# - GMs*(m_jup/m_sun)*(y2-y)/r2**3) #add jupiter mass dependence on asteroid by uncommenting

#position values
# x, y is asteroid
# x2, y2 is jupiter
        for i in range(n.size):
            x[i] = x[i] + vx[i]*h
            y[i] = y[i] + vy[i]*h
        x2 = x2 + vx2*h
        y2 = y2 + vy2*h
   #graph
        xarr = np.append(xarr,x)
        yarr = np.append(yarr,y)
        xarr2 = np.append(xarr2,x2)
        yarr2 = np.append(yarr2,y2)
        count += 1

        t = t+h
        if (count %200 ==0):
            print(f'{count}') #prints count every loop, use to test if it is running and makes it easy to visualize where you are in the code. It will be done at 10,000
    xarr=xarr.reshape([n.size,count+1])
    yarr=yarr.reshape([n.size,count+1])

    r_initial = np.sqrt((xarr[:,0])**2 + (yarr[:,0])**2)
    r_final =np.sqrt((xarr[:,-1])**2 + (yarr[:,-1])**2)
    f =open("asteroid_datafile.txt","a")
    for i in range (n.size):
        if abs(r_final[i] - r_initial[i]) < 1e-1:
            #print(semimajor_a[i])
            f.write(str(semimajor_a[i]) + "\n")
    f.close()

    #now plot
    fig, ax = plt.subplots()

              #plot of jupiter and final positions of asteroids
    ax.plot (xarr[:,-1], yarr[:,-1], 'b', marker='.', linewidth=0)
    ax.plot (xarr2, yarr2, 'r', marker='o')
    ax.set(title='asteroid belt', xlabel='x axis', ylabel='y axis')
    ax.grid()
    fig.savefig('asteroids_final.png') # uncomment to save plot

    plt.show()
#Asteroids() #uncomment to run code and add data to histogram

# def histogram(): #time to plot the histogram 
#     import matplotlib.pyplot as plt
#     import numpy as np
#     import matplotlib

#     f= open('asteroid_datafile.txt','r')
#     a=np.array([])
#     for line in f:
#         s=line.split()
#         a = np.append(a,float(s[0]))
#     f.close()

#     bins = np.linspace(2, 3.6, 500)

#     fig, ax = plt.subplots()
#     plt.hist(a, bins, histtype='bar', color='b')
#     plt.xlabel('AU')
#     plt.ylabel('Number of Asteroids')
#     plt.title('Kirkwood Gaps')
#     fig.savefig('histograms_final.png')
#     plt.show()
# histogram()
