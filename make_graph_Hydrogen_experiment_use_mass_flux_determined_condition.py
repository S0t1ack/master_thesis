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

import numpy as np

topright = False
minorticks = True

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 15
	# axes
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.labelpad'] = 15
	# legend
plt.rcParams['legend.facecolor'] = 'w'
plt.rcParams['legend.edgecolor'] = 'k'
plt.rcParams['legend.framealpha'] = 1.
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.shadow'] = False
	# ticks
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.top'] = topright
plt.rcParams['ytick.right'] = topright
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['xtick.major.width'] = 1.
plt.rcParams['ytick.major.width'] = 1.
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.minor.size'] = 4
if minorticks is True :
	plt.rcParams['xtick.minor.visible'] = True
	plt.rcParams['ytick.minor.visible'] = True
	plt.rcParams['xtick.minor.width'] = 1.0
	plt.rcParams['ytick.minor.width'] = 1.0
	# savefig
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.facecolor'] = 'w'
plt.rcParams['savefig.edgecolor'] = 'w'


def calculate_Re_number(mass_flux, P, subcooled_temperature, fluid, d):
   
    T_main = CP.PropsSI('T', 'P',P,'Q',0,fluid) - subcooled_temperature # 沸騰する温度からサブクール度を引いている
    fluid_rho_liquid = CP.PropsSI('Dmass', 'T|liquid',T_main,'P',P,fluid) # 液相の密度      
    fluid_viscocity = CP.PropsSI('viscosity', 'T', T_main,'P',P,fluid) # 粘度[Pa・s] 
    
    A=100e-6
    velocity = mass_flux*np.pi/4*(15e-3)**2*1/fluid_rho_liquid/A
    
    Re = velocity*d *fluid_rho_liquid/fluid_viscocity
    return Re

def calculate_cavitation_number(mass_flux, P_main, subcooled_temperature, fluid):
    T_main = CP.PropsSI('T', 'P',P_main,'Q',0,fluid) - subcooled_temperature # 沸騰する温度からサブクール度を引いている
    fluid_rho_liquid = CP.PropsSI('Dmass', 'T|liquid',T_main,'P',P,fluid) # 液相の密度   
    fluid_saturated_pressure = CP.PropsSI('P', 'T',T_main,'Q',0,fluid) # 温度時点での飽和蒸気圧
    
    A=100e-6
    velocity = mass_flux*np.pi/4*(15e-3)**2*1/fluid_rho_liquid/A
    
    cavitation_number = (P_main - fluid_saturated_pressure)/(1/2*fluid_rho_liquid*velocity**2)
    return cavitation_number

def return_sigma_and_Re_by_variable(variable,variable_value,subcooled_temperature,fluid,d,mass_flux,P):
    if variable == "Pressure [Pa]":
        sigma = calculate_cavitation_number(mass_flux, variable_value, subcooled_temperature, fluid)
        Re = calculate_Re_number(mass_flux, variable_value, subcooled_temperature, fluid, d)
    elif variable == "Mass Flux [kg/m^2s]":
        sigma = calculate_cavitation_number(variable_value, P, subcooled_temperature, fluid)
        Re = calculate_Re_number(variable_value, P, subcooled_temperature, fluid, d)
    elif variable == "Subcooled Temperature [K]":
        sigma = calculate_cavitation_number(mass_flux, P, variable_value, fluid)
        Re = calculate_Re_number(mass_flux, P, variable_value, fluid, d)        
    
    return sigma,Re


#二分法（方程式の関数項、探索区間の左端、探索区間の右端、キャビテーション数、誤差範囲、最大反復回数）
def bisection(func_f, x_min, x_max, cavitation_number, mass_flux, P,fluid,step, error=1e-10, max_loop=100):
    #初期値を表示
    num_calc = 0  #計算回数
    #print("{:3d}:  {:.15f} <= x <= {:.15f}".format(num_calc, x_min, x_max))
    while True:
        if func_f(mass_flux , P, x_max, fluid)-cavitation_number < 0: 
            x_max+=step 
            if func_f(mass_flux , P, x_max, fluid)-cavitation_number > 0 and func_f(mass_flux , P, x_max-step, fluid)-cavitation_number < 0:
                x_min=x_max-step
                break     
            
        if func_f(mass_flux , P, x_min, fluid)-cavitation_number > 0:
            x_min-=step 
            if func_f(mass_flux , P, x_min+step, fluid)-cavitation_number > 0 and func_f(mass_flux , P, x_min, fluid)-cavitation_number < 0:
                x_max=x_min+step
                break
        
        if func_f(mass_flux , P, x_max, fluid)-cavitation_number > 0 and func_f(mass_flux , P, x_min, fluid)-cavitation_number < 0:
            break
    
    #中間値の定理の条件を満たすか調べる
    if(0 < (func_f(mass_flux , P, x_min, fluid)-cavitation_number)*(func_f(mass_flux , P, x_max, fluid)-cavitation_number)):
        
        print("error: Section definition is invalid (0.0 < func_f(x_min)*func_f(x_max)).")
        #quit()

    #ずっと繰り返す
    while(True):
        #新たな中間値の計算
        x_mid = (x_max +x_min)/2.0

        #探索区間を更新
        if (0 < (func_f(mass_flux , P, x_mid, fluid)-cavitation_number)*(func_f(mass_flux , P, x_max, fluid)-cavitation_number)):  #中間と右端の値が同じの時
            x_max = x_mid  #右端を更新
        else:  #中間と左端の値が同じの時
            x_min = x_mid  #左端を更新

        #結果を表示
        num_calc += 1  #計算回数を数える
        #print("{:2d}:  {:.5f} <= x <= {:.5f}".format(num_calc, x_min, x_max))

        #「誤差範囲が一定値以下」または「計算回数が一定値以上」ならば終了
        if((x_max-x_min <= error) or max_loop <= num_calc):
            break
        
    #最終的に得られた解
    print("T_sub = {:.5f}".format(x_mid))

    return x_mid


fluid = "Hydrogen" 
P=300e3
mass_flux = 130
subcooled_temperature = 0.1
d=5e-3

T= CP.PropsSI('T', 'P',P,'Q',0,fluid) - subcooled_temperature
rho_liquid = CP.PropsSI('Dmass', 'T|liquid',T,'P',P,fluid)

Q = mass_flux*np.pi/4*(15e-3)**2*1/rho_liquid
A=100e-6
velocity = Q/A

Re_list=[]
sigma_list=[]

# --------------------------------------------------------------------


variable = "Mass Flux [kg/m^2s]"
variable_range = [80,130]  

legend = "Subcooled Temperature [K]"
legend_range = [0.1,1.]

# --------------------------------------------------------------------

variable_first = 110

variable_any = 90

# --------------------------------------------------------------------

strMas  = r'$\qquad \qquad \cdot mass flux=$' + rf'${mass_flux:}$'
strVal = r'$\qquad \qquad \cdot variable:$' + rf'${variable:}$'
strFlu = r'$\qquad \qquad \cdot fluid:$' + rf'${fluid:}$'
strRa = r'$\qquad \qquad \cdot range:$' + rf'${variable_range[0]:}$'+' - '+rf'${variable_range[1]:}$'
strPre = r'$\qquad \qquad \cdot Pressure=$' + rf'${P:}$'
strTe = r'$\qquad \qquad \cdot subcooled temperature=$' + rf'${subcooled_temperature:}$'


di = {"Mass Flux [kg/m^2s]":strMas, "Pressure [Pa]":strPre,  "Fluid":strFlu, "Subcooled Temperature [K]":strTe}

legend_di = {"Mass Flux [kg/m^2s]":"m", "Pressure [Pa]":"P",  "Subcooled Temperature [K]":r'$T_{sub}$'}

strlist = [
    strVal,strRa,strFlu, strMas, strPre, strTe
    ]

strlist.remove(di[variable])
strlist.remove(di[legend])

variable_list = np.linspace(variable_range[0],variable_range[1],1000)
legend_list   = np.linspace(legend_range[0],legend_range[1],10)

figure_size=8

Re_figure=plt.figure(figsize=(figure_size,figure_size))
ax_Re=Re_figure.add_subplot(1,1,1)

sigma_figure=plt.figure(figsize=(figure_size,figure_size))
ax_sigma=sigma_figure.add_subplot(1,1,1)


for legend_index in legend_list:
    Re_list.clear()
    sigma_list.clear()
    if legend == "Subcooled Temperature [K]":
        for variable_index in variable_list:
            sigma, Re = return_sigma_and_Re_by_variable(variable, variable_index, legend_index, fluid, d, mass_flux, P)
            Re_list.append(Re)
            sigma_list.append(sigma)
    elif legend == "Pressure [Pa]":
        for variable_index in variable_list:
            sigma, Re = return_sigma_and_Re_by_variable(variable, variable_index, subcooled_temperature, fluid, d, mass_flux, legend_index)
            Re_list.append(Re)
            sigma_list.append(sigma)
    elif legend == "Mass Flux [kg/m^2s]":        
        for variable_index in variable_list:
            sigma, Re = return_sigma_and_Re_by_variable(variable, variable_index, subcooled_temperature, fluid, d, legend_index,P)
            Re_list.append(Re)
            sigma_list.append(sigma)     
            
    ax_Re.plot(variable_list, Re_list,label=legend_di[legend]+"="+format(legend_index, ".1f"))
    ax_sigma.plot(variable_list, sigma_list,label=legend_di[legend]+"="+format(legend_index, ".1f"))

sigma_list.clear()
for legend_index in legend_list:
    sigma_list.append(calculate_cavitation_number(variable_first, P, legend_index, fluid))


step = legend_list[1]-legend_list[0]

subcooled_temperture_list_any=[]
for i in range(len(legend_list)-1):
    subcooled_temperture_list_any.append(bisection(calculate_cavitation_number, legend_list[i], legend_list[i+1] , sigma_list[i],variable_any, P,fluid,step))

print(sigma_list)
print(subcooled_temperture_list_any)

ax_sigma.axvline(x=variable_first,c="k")
ax_sigma.axvline(x=variable_any,c="k")

ax_sigma.plot([variable_first]*(len(subcooled_temperture_list_any)+1),sigma_list,"o",c="b")
ax_sigma.plot([variable_any]*(len(subcooled_temperture_list_any)+1),sigma_list,"o",c="r")

ax_sigma.legend(loc="upper right")
ax_Re.legend(loc="lower right")

ax_Re.set_xlabel(variable)
ax_Re.set_ylabel('Re [-]')

ax_sigma.set_xlabel(variable)
ax_sigma.set_ylabel('$\sigma$ [-]')

for i in range(len(sigma_list)-1):
    ax_sigma.text(variable_any,sigma_list[i],format(subcooled_temperture_list_any[i],".3f"))


rs = 0.0
ts = 0.9
dt = 0.05

for label in strlist :
    ax_Re.text(rs,ts,label,transform=ax_Re.transAxes)
    ax_sigma.text(rs,ts,label,transform=ax_sigma.transAxes)
    ts = ts - dt


    
Re_figure.savefig("Re_vs_"+variable.replace("/","").replace(" ","")+"_legend_"+legend.replace("/","").replace(" ","")+"_"+str(legend_range[0])+"to"+str(legend_range[1])+".png")
sigma_figure.savefig("siguma_vs_"+variable.replace("/","").replace(" ","")+"_legend_"+legend.replace("/","").replace(" ","")+"_"+str(legend_range[0])+"to"+str(legend_range[1])+".png")    
plt.show()


