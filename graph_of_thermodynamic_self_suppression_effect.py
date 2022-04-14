#参照URL
#http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table (公式ドキュメント)
#http://mechlog.jp/materials/201_Thermodynamic_Properties.html (日本語の解説)
import sys
sys.path.append(r'~\Lib\site-packages') # Path to library directory
import CoolProp
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot

import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["mathtext.fontset"] = "stix" 
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 22
plt.rcParams['axes.linewidth'] = 1.3# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線

def calculate_large_sigma(WF,T):
    # ＷＦ:流体の種類
    # T:絶対温度
    
    saturated_pressure=CP.PropsSI('P', 'T',T,'Q',0,WF) #飽和蒸気圧は乾き度が0の圧力である      
    rho_gas    = CP.PropsSI('Dmass', 'T|gas',T,'P',saturated_pressure ,WF) # 飽和蒸気圧下の気相の密度      
    rho_liquid = CP.PropsSI('Dmass', 'T|liquid',T,'P',saturated_pressure,WF) # 飽和蒸気圧下の液相の密度      
    cp_liquid  = CP.PropsSI('Cpmass','T|liquid', T, 'P',saturated_pressure, WF) # 飽和蒸気圧下の液相の比熱容量      
    latent_heat= CP.PropsSI('H', 'P',saturated_pressure,'Q',1,WF)-CP.PropsSI('H', 'P', saturated_pressure,'Q',0,WF)
    # 蒸気潜熱は飽和蒸気圧下の乾き度1と乾き度0でのエンタルピーの差分で定義される
    alpha      = CP.PropsSI('L','T|liquid', T, 'P', saturated_pressure,WF)/(rho_liquid*cp_liquid) # 飽和蒸気圧下の熱拡散率
        
    large_siguma=(latent_heat*rho_gas)**2/(rho_liquid**2*cp_liquid*(T)*alpha**0.5) # 定義通りの計算
    
    return large_siguma 

def generate_large_sigma_list(WF):
    # ＷＦ:流体の種類
    minimum_temperture=CP.PropsSI('T_triple', WF) # 液体の三重点
    maximum_temperture=CP.PropsSI('Tcrit', WF) # 液体の臨界点
    
    #minimum_temperture=CP.PropsSI('T_min', WF) # 液体の三重点
    #maximum_temperture=CP.PropsSI('T_reducing', WF) # 液体の臨界点
    temperture_list = np.arange(minimum_temperture,maximum_temperture,0.5)
    
    large_siguma_list=[]
    
    for T in temperture_list:        
        large_siguma_list.append(calculate_large_sigma(WF, T))
    
    return large_siguma_list,temperture_list

large_siguma_figure=plt.figure(figsize=(12,8))
ax=large_siguma_figure.add_subplot(1,1,1)

hydrogen_large_siguma,tindeh = generate_large_sigma_list('Hydrogen')            
water_large_siguma,tindexw = generate_large_sigma_list('Water')            
nitro_large_siguma,tindexn = generate_large_sigma_list('Nitrogen')            
oxygen_large_siguma,tindexo = generate_large_sigma_list('Oxygen')  
#h_large_siguma,tindexh = generate_large_sigma_list('helium')  

ax.plot(tindeh, hydrogen_large_siguma,label="Hydrogen")          
ax.plot(tindexn, nitro_large_siguma,label="Nitrogen")
ax.plot(tindexo, oxygen_large_siguma,label="Oxygen")
ax.plot(tindexw, water_large_siguma,label="Water")
#ax.plot(tindexh, h_large_siguma,label="Helium")

ax.set_xlabel('Temperature  [K]')
ax.set_ylabel(r'$\Sigma$ [m $s^{-3/2}$]' )
ax.set_yscale('log')
ax.axis([0,700,1,1e7*3])
ax.legend()
ax.grid(which="both")

s=calculate_large_sigma("Water", 303)

s2=calculate_large_sigma("Water", 373)

large_siguma_figure.savefig("large_siguma_vs_T")


plt.show()



