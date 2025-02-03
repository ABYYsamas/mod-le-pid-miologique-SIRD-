import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0,):
    itération = int(T / dt)
    S=[S0]
    I=[I0]
    R=[R0]
    D=[D0]
    
    
    for T in range (1, itération ):
        New_s=S[-1]-Beta*S[-1]*I[-1] *dt
        New_i=I[-1]+Beta*S[-1]*I[-1]-MU*I[-1]-Gamma*I[-1] *dt
        New_r=R[-1]+Gamma*I[-1] *dt
        New_d=D[-1]+MU*I[-1] *dt
        
        S.append(New_s)  
        I.append(New_i)  
        R.append(New_r)  
        D.append(New_d)  
     
    return S, I, R, D

def cost(Beta, Gamma, MU, T, dt, S0, I0, R0, D0):
 data = pd.read_csv("mod-le-pid-miologique-SIRD-\\sird_dataset.csv") 
 SIRD=methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)
 New_mse_s= np.mean((S-data["Susceptibles"]) **2)
 New_mse_i=np.mean((I-data["Infectés"])**2)
 New_mse_r=np.mean((R-data["rétablis"])**2)
 New_mse_d=np.mean((D-data["décés"])**2)
 M_S_E= (New_mse_s,New_mse_d,New_mse_i,New_mse_r)
 return New_mse_s,New_mse_i,New_mse_r,New_mse_d




Beta=0.5
Gamma=0.15
MU=0.15
S0=0.99
I0=0.1 
R0=0
D0=0
T=50
dt=0.1

S,I,R,D=methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)
M_S_E = cost(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)
print("MSE pour les Susceptibles:", M_S_E[0])
print("MSE pour les Infectés:", M_S_E[1])
print("MSE pour les Rétablis:", M_S_E[2])
print("MSE pour les Décès:", M_S_E[3])

temps = np.linspace(0, T, len(S))


plt.figure(figsize=(12, 6))
plt.plot(temps, S, label="susceptible")
plt.plot(temps, I, label="infect")
plt.plot(temps, R, label="Recovered")
plt.plot(temps, D, label="Deceased")
plt.xlabel('Time (days)')
plt.ylabel('Percentage of Population')
plt.show()
