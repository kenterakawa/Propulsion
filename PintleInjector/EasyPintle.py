"""
EasyPintle
ピントルインジェクタの形状から各基本諸元を吐き出すスクリプト
"""



import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import importlib
from scipy import interpolate
from scipy import optimize
import configparser
from matplotlib.font_manager import FontProperties

plt.close('all')


class Pintle:
    def __init__(self, setting_file, reload = False):
        if (reload):  #再読込の際はself.settingの値をそのまま使う
            pass
        else:
            self.setting_file = setting_file
            self.setting = configparser.ConfigParser()
            self.setting.optionxform = str  # 大文字小文字を区別するおまじない
            self.setting.read(setting_file, encoding='utf8')
        setting = self.setting

        self.name = setting.get("全体", "名前")
        #self.is_save_fig = setting.getboolean("全体", "図を保存する？")
        Dist_Dp = setting.getfloat("PintleDim", "Pintle Tip Diameter Dp[mm]")/10**3
        Dist_Ls = setting.getfloat("PintleDim", "Skip Distance Ls[mm]")/10**3
        Dist_Dfo = setting.getfloat("PintleDim", "Fuel channel Outer Diameter Dfo[mm]")/10**3
        Dist_Doi = setting.getfloat("PintleDim", "LOx channel Inner Diameter Doi [mm]")/10**3
        Dist_Lo1 = setting.getfloat("PintleDim", "Primary Oxidizer slot height Lo1[mm]")/10**3
        Dist_deltao1 = setting.getfloat("PintleDim", "Primary Oxidizer slot Circumferential size deltao1[mm]")/10**3
        Dist_Lo2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot height Lo2[mm]")/10**3
        Dist_deltao2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot Circumferential size deltao2[mm]")/10**3
        CN_N1 = setting.getfloat("PintleDim", "Number of Primary slots N1")
        CN_N2 = setting.getfloat("PintleDim", "Number of Secondary slots N2")
        Dist_Lo12 = setting.getfloat("PintleDim", "Distance between Primary and Secondary Oxidizer slots Lo12[mm]")/10**3
        rho_o = setting.getfloat("OxidizerProp", "LOx Density[kg/m3]")
        rho_go = setting.getfloat("OxidizerProp", "LOx Two Phase Density[kg/m3]")
        mdot_o = setting.getfloat("OxidizerProp", "LOx mass flow rate [kg/s]")
        rho_f = setting.getfloat("FuelProp", "Fuel Density[kg/m3]")
        rho_gf = setting.getfloat("FuelProp", "Fuel Two Phase Density[kg/m3]")
        mdot_f = setting.getfloat("FuelProp", "Fuel mass flow rate [kg/s]")
        LOx_Cd = setting.getfloat("InjectorParam", "LOx injector Cd")
        Fuel_Cd = setting.getfloat("InjectorParam", "Fuel injector Cd")
        ORD_FMR_a = setting.getfloat("InjectorParam", "Velocity a for FMR [m/s]")

        #Injector geometry area
        self.Area_LOx1 = Dist_Lo1*Dist_deltao1*CN_N1
        self.Area_LOx2 = Dist_Lo2*Dist_deltao2*CN_N2
        self.Area_LOx = self.Area_LOx1+self.Area_LOx2
        self.Area_Fuel = (Dist_Dfo**2-Dist_Dp**2)*np.pi/4
        self.BLF1 = CN_N1*Dist_deltao1/np.pi/Dist_Dp
        self.BLF2 = CN_N2*Dist_deltao2/np.pi/Dist_Dp
        self.LOx_PA_ratio = self.Area_LOx1/(self.Area_LOx1+self.Area_LOx2)

        #flow velocity
        self.LOx_vo1 = mdot_o*self.LOx_PA_ratio/rho_o/self.Area_LOx1
        self.LOx_vo2 = mdot_o*(1-self.LOx_PA_ratio)/rho_o/self.Area_LOx2
        self.LOx_vo = self.LOx_vo1
        self.Fuel_vf = mdot_f/rho_f/self.Area_Fuel
        self.Mom_Ff = rho_f*self.Fuel_vf**2*self.Area_Fuel
        self.Mom_Fo1 = rho_o*self.LOx_vo1**2*self.Area_LOx1
        self.Mom_Fo2 = rho_o*self.LOx_vo2**2*self.Area_LOx2
        self.Mom_Fo = self.Mom_Fo1+self.Mom_Fo2
        self.TMR = self.Mom_Ff/self.Mom_Fo

        #Other Parameters
        self.Skip_Dist_D = Dist_Ls/Dist_Dp
        self.Skip_Dist_V = Dist_Ls/self.Fuel_vf
        self.ATM_Cone = np.arctan(self.TMR**0.5)*180/np.pi

        #Index a for mixing paramter
        #self.mppar_a = (self.Fuel_vf*rho_f/rho_gf*mdot_f+self.LOx_vo*rho_o/rho_go*mdot_o)/(mdot_f+mdot_o)
        
        #Manual Setting (a=184~242 m/s)
        self.mppar_a = ORD_FMR_a
        
        #mixing paramter, final momemtum ratio
        self.mppar_C1 = Dist_Lo1/self.Fuel_vf*self.mppar_a/np.pi/Dist_Dp*self.BLF1/(1-self.BLF1)
        self.mppar_C2 = Dist_Lo2/self.Fuel_vf*self.mppar_a/np.pi/Dist_Dp*self.BLF2/(1-self.BLF2)
        self.mp1 = (rho_f*self.Fuel_vf**2*(Dist_Dfo-Dist_Dp)/2*(Dist_deltao1+2*self.mppar_C1*Dist_Lo1))/(rho_o*self.LOx_vo1**2*Dist_deltao1*Dist_Lo1)
        self.mp2 = (rho_f*self.Fuel_vf**2*(Dist_Dfo-Dist_Dp)/2*(Dist_deltao2+2*self.mppar_C2*Dist_Lo2))/(rho_o*self.LOx_vo1**2*Dist_deltao2*Dist_Lo2)

        #Pressure loss delta p
        self.deltap_o =  mdot_o**2/(2*rho_o*self.Area_LOx**2*LOx_Cd**2)/10**6       
        self.deltap_f = mdot_f**2/(2*rho_f*self.Area_Fuel**2*Fuel_Cd**2)/10**6
        

    def display(self):
    	 print("")
    	 print("TMR (Total Momentum Ratio) :\t\t%.2f " % (self.TMR))    	 
    	 print("LOx Primary slot flow ratio:\t\t%.2f " % (self.LOx_PA_ratio))
    	 print("LOx Outlet Velocity vo1 [m/s] :\t\t%.1f " % (self.LOx_vo1))
    	 print("Fuel Outlet Velocity vf [m/s] :\t\t%.1f " % (self.Fuel_vf))    	 
    	 print("")
    	 print("Non-Dimensional skip distance Ls/Dp :\t%.2f " % (self.Skip_Dist_D))
    	 print("Normalized skip distance Ls/vf[s] :\t%.5f " % (self.Skip_Dist_V))
    	 print("Theoretical Atomizing Cone angle[deg]:\t%.2f " % (self.ATM_Cone))
    	 print("")
    	 print("Primary slot Blockage Factor :\t\t%.2f " % (self.BLF1))
    	 print("Secondary slot Blockage Factor :\t%.2f " % (self.BLF2))
    	 print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f))
    	 print("Oxidizer Injector delta p [MPa]:\t%.2f " % (self.deltap_o))
    	 print("")
    	 #print("Primary Mixing Parameter coefficient a:\t\t%.2f " % (self.mppar_a))
    	 print("Primary Mixing Parameter :\t\t%.2f " % (self.mp1)) 
    	 print("Secondary Mixing Parameter :\t\t%.2f " % (self.mp2))
    	 print("Primary Mixing Parameter coefficient C1:\t\t%.2f " % (self.mppar_C1))

    def print(self):
    	with open("PintleParams.out","w") as output:
    	 print("TMR (Total Momentum Ratio) :\t\t%.2f " % (self.TMR),file=output)
    	 print("LOx Primary slot flow ratio:\t\t%.2f " % (self.LOx_PA_ratio),file=output)
    	 print("LOx Outlet Velocity vo1 [m/s] :\t\t%.1f " % (self.LOx_vo1),file=output)
    	 print("Fuel Outlet Velocity vf [m/s] :\t\t%.1f " % (self.Fuel_vf),file=output)
    	 print("Non-Dimensional skip distance Ls/Dp :\t%.2f " % (self.Skip_Dist_D),file=output)
    	 print("Normalized skip distance Ls/vf[s] :\t%.5f " % (self.Skip_Dist_V),file=output)
    	 print("Primary slot Blockage Factor :\t\t%.2f " % (self.BLF1),file=output)
    	 print("Secondary slot Blockage Factor :\t%.2f " % (self.BLF2),file=output)
    	 print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f),file=output)
    	 print("Oxidizer Injector delta p [MPa]:\t%.2f " % (self.deltap_o),file=output)
    	 print("Primary Mixing Parameter :\t\t%.2f " % (self.mp1),file=output) 
    	 print("Secondary Mixing Parameter :\t\t%.2f " % (self.mp2),file=output) 
    	 #print("Primary Mixing Parameter coefficient a:\t\t%.2f " % (self.mppar_a),file=output) 




if __name__ == '__main__':
    if len(sys.argv) == 1:
        setting_file = 'setting.ini'
    else:
        setting_file = sys.argv[1]
        assert os.path.exists(setting_file), "ファイルが存在しません"
    plt.close("all")
    plt.ion()
    pout = Pintle(setting_file)
    pout.display()
    pout.print()

 
