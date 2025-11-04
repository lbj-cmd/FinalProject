#Cassandra Bodin
#Phys305
#Final Project: Asteroid Version

#run by "python final_project_asteroids_bodin"

#Goal of Project:
   #Recreate the threebody problem we did for homework, for Jupiter the Sun and asteroids.
   #Instead of just 4 asteroids though I want to add a large amount of ateroids, maybe 100-1000, all with different semi major axes lengths
   #Plot and animate each in different colors
   #Add code to take the data for the asteroids and use it to plot a histogram showing the Kirkwood Gaps

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

def histogram(): #time to plot the histogram 
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib

    f= open('asteroid_datafile.txt','r')
    a=np.array([])
    for line in f:
        s=line.split()
        a = np.append(a,float(s[0]))
    f.close()

    bins = np.linspace(2, 3.6, 500)

    fig, ax = plt.subplots()
    plt.hist(a, bins, histtype='bar', color='b')
    plt.xlabel('AU')
    plt.ylabel('Number of Asteroids')
    plt.title('Kirkwood Gaps')
    fig.savefig('histograms_final.png')
    plt.show()
histogram()
