# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 22:02:08 2020

@author: Zihao Li
"""

from matplotlib.pyplot import*
from numpy import*
import mpl_toolkits.axisartist as axisartist
import imageio
import glob

def plotwave(m,a,V0,E):
    x0 = 30e-10
    hbar = 1.055e-34
    x1 = linspace(-x0,0,500)
    x2 = linspace(0,a,100)
    x3 = linspace(a,x0,500)
    
    k1 = (2*mu*E/hbar**2)**0.5
    k3 = (2*mu*(V0-E)/hbar**2)**0.5
    k2 = 1j*k3
    
    A1 = 1
    
    m = e**(1j*k2*a)
    n = e**(-1j*k2*a)
    q = e**(1j*k1*a)

    C = -4*A1*k1*k2*m*n/(k1**2*m-2*k1*k2*m+k2**2*m-k1**2*n-2*k1*k2*n-k2**2*n)/q
    A2 = A1*(k1**2*m-k2**2*m-k1**2*n+k2**2*n)/(k1**2*m-2*k1*k2*m+k2**2*m-k1**2*n-2*k1*k2*n-k2**2*n)

    
    B1 = (2*A1*k1*(k1 + k2)*n)/(-k1**2*m + 2*k1*k2*m - k2**2*m + k1**2*n + 2*k1*k2*n + k2**2*n)
    B2 = (2*A1*k1*(k2 - k1)*m)/(-k1**2*m + 2*k1*k2*m - k2**2*m + k1**2*n + 2*k1*k2*n + k2**2*n)
    
    y1 = A1*e**(1j*k1*x1)+A2*e**(-1j*k1*x1)
    y2 = B1*e**(1j*k2*x2)+B2*e**(-1j*k2*x2)
    y3 = C*e**(1j*k1*x3)
    

    fig = figure(figsize = (6,4.5))
    
    ax1 = axisartist.Subplot(fig, 111)
    
    fig.add_axes(ax1)
    ax1.plot(x1,y1,'b')
    ax1.plot(x2,y2,'b')
    ax1.plot(x3,y3,'b')
    ax1.set_ylim(-2,2)
    ax1.set_xlim(-1.05*x0,1.05*x0)
    
    
    ax1.axis["x"] = ax1.new_floating_axis(0,0)
    ax1.axis["x"].set_axisline_style("-|>", size = 2)
    ax1.axis["y"] = ax1.new_floating_axis(1,0)
    ax1.axis["y"].set_axisline_style("-|>", size = 2)
    ax1.axis["x"].set_label('x/(nm)')
    ax1.set_yticks([])
    ax1.axis["top"].set_visible(False)
    ax1.axis["right"].set_visible(False)
    ax1.axis["bottom"].set_visible(False)
    ax1.axis["left"].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    
    ax1.text(1.5e-9, 1.2, r'$\mu$=%skg'%format(mu,'.2e')+
             '\n'+r'$V_0$=%.4f$eV$'%(V0/1.6021766208/10**(-19))+
             '\n'+r'E=%.4f$eV$'%(E/1.6021766208/10**(-19))+
             '\n'+r'a=%.4f$\AA$'%(a*1e10)+'\n'+
             r'D=%.4f'%abs(C)**2+'\n'+
             r'R=%.4f'%abs(A2)**2)

    ax1.vlines(a,0,0.5*V0/1.6021766208/10**(-19),'black','--')
    ax1.vlines(0,0,0.5*V0/1.6021766208/10**(-19),'black','--')
    ax1.hlines(0.5*V0/1.6021766208/10**(-19),0,a,'black','--')

# a = 1e-10
# V0 = 2*1.6021766208*10**(-19)
# E = 1.6021766208*10**(-19)
# mu = 9.10956*10**(-31)
# plotwave(mu,a,V0,E)
    
n = 30
i = 1

# 以下4个循环分别控制势垒宽度a,粒子质量mu,势垒能量V0,粒子能量E，可根据需求自行调节
# 输出动图需要分别运行以下四个循环
for a in linspace(0.01e-10,15e-10,n):
    
    V0 = 2*1.6021766208*10**(-19)
    E = 1*1.6021766208*10**(-19)
    mu = 9.10956*10**(-31)
    
    plotwave(mu,a,V0,E)
    savefig('%d.png'%(i))
    close()
    i = i+1


# for mu in linspace(9.10956*10**(-31),50*9.10956*10**(-31),n):
    
#     V0 = 2*1.6021766208*10**(-19)
#     E = 1.6021766208*10**(-19)
#     a = 2e-10
    
#     plotwave(mu,a,V0,E)
#     savefig('%d.png'%(i))
#     close()
#     i = i+1

# for V0 in linspace(1.5*1.6021766208*10**(-19),10*1.6021766208*10**(-19),n):
    
#     E = 1.6021766208*10**(-19)
#     a = 2e-10
#     mu = 9.10956*10**(-31)
    
#     plotwave(mu,a,V0,E)
#     savefig('%d.png'%(i))
#     close()
#     i = i+1


# for E in linspace(0.1*1.6021766208*10**(-19),1.99*1.6021766208*10**(-19),n):
    
#     V0 = 2*1.6021766208*10**(-19)
#     a = 2e-10
#     mu = 9.10956*10**(-31)
    
#     plotwave(mu,a,V0,E)
#     savefig('%d.png'%(i))
#     close()
#     i = i+1

def create_gif(image_list, gif_name):  
  
    frames = []  
    for image_name in image_list:  
        frames.append(imageio.imread(image_name))  
    # Save them as gif files  
    imageio.mimsave(gif_name, frames, 'GIF', duration = 0.05)  
  
    return 


def find_all_png():
    buf=[]
    for i in range(int(n)):
        png_filenames = glob.glob('%d.png'%(i))
        for png_file in png_filenames:
            buf.append(png_file)
    return buf


if __name__ == '__main__':#Get gif file
 
    buff = find_all_png()
    create_gif(buff,'result.gif' ) 