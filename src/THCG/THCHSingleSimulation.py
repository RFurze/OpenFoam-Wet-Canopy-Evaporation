import numpy as np
import matplotlib.pyplot as plt

#Rain fall event

total_rain_fall=30 #mm
time=3 #hours
r=total_rain_fall/(time*60)

#inital conditions


#coeffecents

per=1
cp=1
cv=1
L=2260
rho_a=1.293
P_a=101000
r_a=0.2
h=20
gamma=1
beta=P_a/(rho_a*h)
eta=0.622
h_c=(rho_a*cp)/(r_a)
alpha = 0.9*per   #mm per min
c_w = 1.36*per    #mm
a = 3.7     #mm per min
d = 0.0002    #mm per min
S = c_w       #mm per min


#starting values
Rn=5
H_o = 0.3     # %
T_o = 21      # c^o
T_co = 21     # c^o
T_amb = 20    # c^o
v=0        #mm
u=0       #mm

#rain input function
def R(t,r):
    
    if 30<t<(30+time*60):    
        fun = r  
    else:
        fun =0       
    return fun

#function for saturation pressure
def e_s(T):
    return (np.exp(34.494-(4924.99/(T+237.1))))/((T+105)**1.52)#https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml?tab_body=pdf

e_a=H_o*e_s(T_amb)
print(e_a)
e=per*rho_a*eta/(P_a*r_a)
#function for evaporation
def E(e_b,T):
    return e*(e_s(T)-e_b)


#stem flow function
def ben(C,c_w=c_w,a=a):    
    if C>c_w:
        fun = d*np.exp(a*(C-c_w))    
    else:
        fun=0
    
    return fun

def sat(C,S=S):
    if C<S:
        Sat=C/S
    else:
        Sat=1
    return Sat



def lorenz(C,G,e_b,T,Tc,t, R=R,  c_w=c_w, T_amb=T_amb):
    '''
    Given:
       C,G, t: a point of interest in three dimensional space
       a, b, c, d, e, f, g, h: parameters defining the lorenz attractor
    Returns:a
       C_dot,G_dot, t_dot: values of the lorenz attractor's partial
           derivatives at the point C, G, t
    '''

    C_dot = alpha*R(t,r)-1.003*E(e_b,Tc)*sat(C,S)-ben(C,c_w,a)
    
    G_dot = (1-alpha)*R(t,r) +ben(C,c_w,a)
    
    e_b_dot = beta*E(e_b,Tc)*sat(C,S)-gamma*(e_b-e_a)
    
    T_dot = h_c*(Tc-2*T+T_amb)/20*cp 
    
    Tc_dot = (Rn-h_c*(Tc-T)-L*E(e_b,Tc)*sat(C,S))/200*cv
    
    t_dot = 1
    return C_dot,G_dot,e_b_dot,T_dot,Tc_dot, t_dot



dt = 0.01
window=time+3 #hours
num_steps = (window*60)*100


# Need one more for the initial values
#creates arrays

Cs = np.empty((num_steps + 1,))
Gs = np.empty((num_steps + 1,))
ts = np.empty((num_steps + 1,))
dC = np.empty((num_steps + 1,))
dy = np.empty((num_steps + 1,))
dr = np.empty((num_steps + 1,))
rs = np.empty((num_steps + 1,))
e_b = np.empty((num_steps + 1,))
T = np.empty((num_steps + 1,))
Tc = np.empty((num_steps + 1,))

  

# Set initial values

rs[0],dr[0],Cs[0],Gs[0],ts[0],dC[0],dy[0],e_b[0],T[0],Tc[0] = (0,0,v,u,0,0,0, H_o*e_s(T_co),T_o,T_co)


# Step through "time", calculating the partial derivatives at the current point
# and using them to estimate the neCt point

for i in range(num_steps):
    C_dot, G_dot, e_b_dot, T_dot, Tc_dot, t_dot = lorenz(Cs[i], Gs[i], e_b[i], T[i],Tc[i],ts[i])
    dr[i + 1] = R((i+1)*dt,r)
    Cs[i + 1] = Cs[i] + (C_dot * dt)
    Gs[i + 1] = Gs[i] + (G_dot*dt)
    ts[i + 1] = ts[i] + (t_dot * dt)
    
    dC[i + 1] = C_dot
    dy[i + 1] = G_dot
    rs[i + 1] = rs[i]+dr[i+1]*dt
    
    e_b[i + 1] = e_b[i]+e_b_dot*dt
    T[i + 1] = T[i]+T_dot*dt
    Tc[i + 1] = Tc[i]+Tc_dot*dt
    
 

# Plot solutions
ts=ts/60

evap = np.empty((num_steps+1 ,))
drain = np.empty((num_steps+1 ,))
H = np.empty((num_steps+1 ,))
for i in range(num_steps):
    evap[i]=E(e_b[i],Tc[i])*sat(Cs[i],S)
    drain[i]=ben(Cs[i],c_w,a)
    H[i]=e_b[i]/e_s(Tc[i])
evap[num_steps]=evap[num_steps-1]
drain[num_steps]=drain[num_steps-1]
H[num_steps]=H[num_steps-1]

#plt.figure()
f, axes = plt.subplots(3, 1)

#axes[0].subplot(3,1,1)
axes[0].plot(ts, rs-(Gs+Cs)  ,'skyblue', label='Evapiration')
axes[0].plot(ts, rs , 'k-',linestyle='dashed', label='Rain Fall')
axes[0].plot(ts, Cs  , 'g-', label='Canopy')
axes[0].plot(ts,Gs  , 'b-', label='Ground')

axes[0].grid()
axes[0].legend(loc='best')
axes[0].set_ylabel('Water present (mm)')
axes[0].set_xlim([-0.25, ts[i]+1.25])



#plt.subplot(3,1,2)
axes[1].plot(ts,evap,'skyblue', label='Evapiration')
axes[1].plot(ts,dr,'k-',linestyle='dashed',label='Rain Fall')
axes[1].plot(ts, dC  , 'g-', label='Canopy')
axes[1].plot(ts, dy  , 'b-', label='Ground')
axes[1].grid()
axes[1].set_ylabel('Rate of change (mm/min)')
axes[1].set_xlim([-0.25, ts[i]+1.25])
axes[1].legend(loc='best')


#plt.subplot(3,1,3)
axes[2].plot([-1,-1],[10,25],'k-',label='Humidity',linestyle='dashed')
axes[2].plot(ts,T,'r-',label='Temp air')
axes[2].plot(ts,Tc,'c-',label='Temp leaves')
axes[2].plot([0,ts[num_steps]],[T_amb,T_amb],'gray',label='Amb Temp',linestyle='dashed')
axes[2].set_xlim([-0.25, ts[i]+1.25])
axes[2].grid()
axes[0].plot(ts, rs-(Gs+Cs)  ,'skyblue', label='Evapiration')
axes[2].set_xlabel('Time (hrs)')
axes[2].legend(loc='best')
axes[2].set_ylabel('Temprature (C^o)')
ax2=axes[2].twinx()
ax2.set_ylabel('Humidity (%)')
ax2.plot(ts, 100*H  , 'k-',linestyle='dashed')
ax2.plot([-1,-1],[0,100],'k-',label='Humidity',linestyle='dashed')

plt.show()


#quantifies temporal effects

print('######################################################################')
print('total rain fall =',rs[num_steps],'mm')
print('total water in ground =',Gs[num_steps],'mm')
print('total water evaporated =',rs[num_steps]-Gs[num_steps],'mm')
sim_loss=100*(rs[num_steps]-Gs[num_steps])/rs[num_steps]
print('percentage of water evaporated =',sim_loss,'%')
print('######################################################################')
est_loss=100*np.max(1.003*evap)/(a*r)
print('evaporation/rate of raine fall prediction=',est_loss,'%')
print('evaporation during drip aditonal contrabution =',sim_loss-est_loss,'% ') 
print('contrabution to over all effect =',100*(sim_loss-est_loss)/sim_loss,'%')
print('######################################################################')
print('stedy state evaporation rate=', evap[14500],'mms^-1')
print('total precipitation=',time*r*60,'mm')
print('E/R ratio=',evap[14500]/r)
