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
        self.PTDia = setting.getfloat("PintleDim", "Pintle Tip Diameter Dp[mm]")
        self.Dist_Ls = setting.getfloat("PintleDim", "Skip Distance Ls[mm]")
        self.Dist_Dfo = setting.getfloat("PintleDim", "Fuel channel Outer Diameter Dfo[mm]")
        self.Dist_Doi = setting.getfloat("PintleDim", "LOx channel Inner Diameter Doi [mm]")
        self.Dist_Lo1 = setting.getfloat("PintleDim", "Primary Oxidizer slot height Lo1[mm]")
        self.Dist_deltao1 = setting.getfloat("PintleDim", "Primary Oxidizer slot Circumferential size deltao1[mm]")
        self.Dist_Lo2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot height Lo2[mm]")
        self.Dist_deltao2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot Circumferential size deltao2[mm]")
        self.SP_N1 = setting.getfloat("PintleDim", "Number of Primary slots N1")
        self.SP_N2 = setting.getfloat("PintleDim", "Number of Secondary slots N2")
        self.Dist_Lo12 = setting.getfloat("PintleDim", "Distance between Primary and Secondary Oxidizer slots Lo12[mm]")
        self.rho_o = setting.getfloat("OxidizerProp", "LOx Density[kg/m3]")
        self.mdot_o = setting.getfloat("OxidizerProp", "LOx mass flow rate [kg/s]")
        self.rho_f = setting.getfloat("FuelProp", "Fuel Density[kg/m3]")
        self.mdot_f = setting.getfloat("FuelProp", "Fuel mass flow rate [kg/s]")
        self.LOx_Cv = setting.getfloat("InjectorParam", "LOx injector Cv")
        self.Fuel_Cv = setting.getfloat("InjectorParam", "Fuel injector Cv")

        #LOx side injector area
        self.Area_LOx = 1



    def display(self):
    	 print("TMR (Total Momentum Ratio) :\t\t%.1f " % (self.Area_LOx))

    def print(self):
    	with open("PintleParams.out","w") as output:
    		print("TMR (Total Momentum Ratio) :\t\t%.1f " % (self.Area_LOx),file=output)  



"""
CL出力
    def print(self):
        with open("tankout_S1.out","w") as output:
            print("タンク鏡重量:\t%.1f [kg]" %(self.hoge),file=output) #add
       
"""

if __name__ == '__main__':
    if len(sys.argv) == 1:
        setting_file = 'setting.ini'
    else:
        setting_file = sys.argv[1]
        assert os.path.exists(setting_file), "ファイルが存在しません"
    plt.close("all")
    plt.ion()
    pintout = Pintle(setting_file)
    pintout.display()
    #pintout.print()
 
