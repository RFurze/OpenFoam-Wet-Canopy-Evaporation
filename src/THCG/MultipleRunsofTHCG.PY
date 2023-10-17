import numpy as np
import matplotlib.pyplot as plt

#Rain fall event

#total_rain_fall=30 #mm



#inital conditions




def sim(A):
    Rt=0.1
    r=A[0]
    X=[]
    Y=[]
    rs=[1]
    num_steps=0
    c_w = A[1]    #mm
    alpha =  A[2] #mm per min
    gamma=A[3]
    a = A[4]     #mm per min
    H_o = A[5]     # %
    while Rt<121:
        if Rt>10:
            mark=1
        else:
            mark=1
        time=int(Rt/r)
        #print(time)
        #print(Rt)
        per=1
        cp=1
        cv=1
        L=2260
        rho_a=1.293
        P_a=101000
        r_a=0.2
        h=20
        
        beta=P_a/(rho_a*h)
        eta=0.622
        h_c=(rho_a*cp)/(r_a)
        d = 0.0002    #mm per min
        S = c_w       #mm per min
        
        
        #starting values
        Rn=0
        
        T_o = 21      # c^o
        T_co = 21     # c^o
        T_amb = 20    # c^o
        v=0        #mm
        u=0       #mm
        
        #rain input function
        def R(t,r):
            
            if 180<t<(180+time):    
                fun = r  
            else:
                fun =0       
            return fun
        
        #function for saturation pressure
        def e_s(T):
            return (np.exp(34.494-(4924.99/(T+237.1))))/((T+105)**1.52)#https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml?tab_body=pdf
        
        e_a=H_o*e_s(T_amb)
        e=per*rho_a*eta/(P_a*r_a)
        #function for evaporation
        def E(e_b,T,e):
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
        
            C_dot = alpha*R(t,r)-1.003*E(e_b,Tc,e)*sat(C,S)-ben(C,c_w,a)
            
            G_dot = (1-alpha)*R(t,r) +ben(C,c_w,a)
            
            e_b_dot = beta*E(e_b,Tc,e)*sat(C,S)-gamma*(e_b-e_a)
            
            T_dot = h_c*(Tc-2*T+T_amb)/20*cp 
            
            Tc_dot = (Rn-h_c*(Tc-T)-L*E(e_b,Tc,e)*sat(C,S))/200*cv
            
            t_dot = 1
            return C_dot,G_dot,e_b_dot,T_dot,Tc_dot, t_dot
        
        
        
        dt = 1/mark
        window=time+10000 #hours
        num_steps = int((window)*mark)
        
        
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
    
    
        sim_loss=100*(rs[num_steps]-Gs[num_steps])/rs[num_steps]
        Y.append(sim_loss)
        X.append(rs[num_steps])

        if Rt<10:
            Rt=Rt+0.2
        else:
            Rt=Rt+1
 
    loss=round(E(e_b[time],Tc[time],e)*sat(Cs[time],S)/r,2)
    jeez=round(c_w,2)
    oss=round(alpha,1)
    
    plt.plot(X,Y,label=['E/P=',loss,' c_w=',jeez,'mm',' Alpha=',oss])
    #print(rs[num_steps])
############################################################################################
i=1
while i<11:
    l=i/10
    A=[0.08,1.3*(1-l)+0.3,1.1-l,0.1,0.4,l]
    sim(A)
    i=i+1
        


'''
r=A[0]=
c_w = A[1]    
alpha =  A[2] 
gamma=A[3]
a = A[4]
H_o = A[5]
'''


plt.grid()
plt.ylim([0,100.3])
plt.xlim([0,120])
plt.title('THCG Model Plot')
plt.ylabel('Ewc (%of Pg)')
plt.xlabel('Pg (mm)')
plt.legend(loc='best')



plt.show()
