import sys
import numpy as np
import cmath
from math import exp, pi, atan
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
from scipy.sparse import bsr_array
import ROOT
plt.rc('text', usetex=False)
from histos import Hist2D, Hist3D, plot_th1, errorfill, plot_tf1

#l, m = int(sys.argv[1]), int(sys.argv[2])

def Y_theta(eta):
    #return 2*np.arctan(np.exp(-eta))
    ## the Spherical harmonic functions has the Legendre component P_l(cos(theta)). In high-energy collisions we use
    ## eta to define the polar angle. It is better introduced in the Legendre functions as P_lm(eta/eta_max) according to 
    ## studies such as arXiv:1506.03496 and arXiv:1909.03979
    ## In SyPy the spherical harmonics are called as sph_harm(m, l, phi, theta). That means, the Legendre function is called as
    ## P_lm(cos(theta)). Lets define theta as such that cos(theta)=eta/eta_max
    eta_max = 6
    return np.arccos(eta/eta_max)

def SH(x, pars):
    #par[0] = a00
    #par[1] = a1-1
    #par[2] = a10
    #par[3] = a11
    #par[4] = a2-2
    #par[5] = a2-1
    #par[6] = a20
    #par[7] = a21
    #par[8] = a22
    #par[9] = a3-3
    #par[10] = a3-2
    #par[11] = a3-1
    #par[12] = a30
    #par[13] = a31
    #par[14] = a32
    #par[15] = a33
    phi = x[0]+pi
    theta = Y_theta(x[1])
    total = 0.0
    ipar = 0
    for l in range(0,4):
        for m in range(-l,l+1):
            Y = sph_harm(abs(m), l, phi, theta) ## keep this sequence. In ScyPy theta is unconventionally the azimuthal angle and phi is the polar angle
            if m < 0: Y = np.sqrt(2) * (-1)**m * Y.imag
            elif m > 0: Y = np.sqrt(2) * (-1)**m * Y.real
            total += pars[ipar]*Y
            ipar += 1
    #total *= complex(total.real, -total.imag)
    return total.real

sh = ROOT.TF2('sh', SH, -pi, pi, -6, 6, 16)

# Grids of polar and azimuthal angles
etabins = np.linspace(-6, 6, 50)
phi = np.linspace(-pi, pi, 50)
# Create a 2-D meshgrid of (theta, phi) angles.
eta, phi = np.meshgrid(etabins, phi)

    
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(Y_theta(eta)) * np.sin(phi),
                np.sin(Y_theta(eta)) * np.cos(phi),
                np.cos(Y_theta(eta))])

def plot_histo(ax,h) :
    scale = h.Integral()
    h.Scale(1/scale)
    data = []
    for x in range(1,h.GetNbinsX()+1):
        datay = []
        for y in range(1,h.GetNbinsY()+1):
            datay.append(h.GetBinContent(x,y))
        data.append(datay)
    data = np.array(data)
    
    theta_ = Y_theta(etabins)
    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    cmap.set_clim(0, data.max())
    ax.plot_surface(phi, etabins, data,
                    facecolors=cmap.to_rgba(data),
                    rstride=2, cstride=2)
    ax.set_xlabel(r'$\Delta\phi$', fontsize='xx-large')
    ax.set_ylabel(r'$\Delta\eta$', fontsize='xx-large')


def plot_histo_xyz(ax,h) :
    scale = h.Integral()
    h.Scale(1/scale)
    data = []
    for x in range(1,h.GetNbinsX()+1):
        datay = []
        for y in range(1,h.GetNbinsY()+1):
            datay.append(h.GetBinContent(x,y))
        data.append(datay)
    data = np.array(data)
    
    theta_ = Y_theta(etabins)
    Yx, Yy, Yz = np.abs(data) * xyz
    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    cmap.set_clim(data.min(), data.max())

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Yz),
                    rstride=2, cstride=2)    
    # Draw a set of x, y, z axes for reference.
    ax_lim = data.max()
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.
    ax.set_title(r'data')
    ax_lim = data.max()
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')

def plot_Y(ax, el, m):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""

    # NB In SciPy's sph_harm function the azimuthal coordinate, theta,
    # comes before the polar coordinate, phi.
    theta_ = Y_theta(etabins) #np.linspace(0, pi, 50)
    phi = np.linspace(-pi, pi, 50)
    theta_, phi = np.meshgrid(theta_, phi)
    Y = sph_harm(abs(m), el, phi, theta_)

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    #Yx, Yy, Yz = np.abs(Y) * xyz

    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    cmap.set_clim(0, Y.real.max())

    ax.plot_surface(phi, theta_, Y.real,
                    facecolors=cmap.to_rgba(Y.real),
                    rstride=2, cstride=2)
    
    ax.set_title(r'$Y_{{{},{}}}$'.format(el, m), fontsize='xx-large') 
    ax.axis('off')

def plot_sh(ax, sh):
    theta_ = Y_theta(etabins) #np.linspace(0, pi, 50)
    phi_ = np.linspace(-pi, pi, 50)
    theta_, phi_ = np.meshgrid(theta_, phi_)
    total =  sh.GetParameter(0)*sph_harm(0, 0, phi_, theta_)
    ipar = 1
    for l in range(1,4):
        for m in range(-l,l+1):
            Y = sph_harm(abs(m), l, phi_, theta_) ## keep this sequence. In ScyPy theta is unconventionally the azimuthal angle and phi is the polar angle
            if m < 0: Y = np.sqrt(2) * (-1)**m * Y.imag
            elif m > 0: Y = np.sqrt(2) * (-1)**m * Y.real
            total += sh.GetParameter(ipar)*Y
            ipar += 1

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    #Yx, Yy, Yz = np.abs(Y) * xyz

    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    cmap.set_clim(0, total.real.max())

    ax.plot_surface(phi_, theta_, total.real,
                    facecolors=cmap.to_rgba(total.real),
                    rstride=2, cstride=2)
            

def plot_Y_xyz(ax, el, m):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""

    # NB In SciPy's sph_harm function the azimuthal coordinate, theta,
    # comes before the polar coordinate, phi.
    theta_ = Y_theta(etabins) #np.linspace(0, pi, 50)
    phi_ = np.linspace(-pi, pi, 50)
    theta_, phi_ = np.meshgrid(theta_, phi_)
    Y = sph_harm(abs(m), el, phi, theta_)

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    Yx, Yy, Yz = np.abs(Y) * xyz

    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('PRGn'))
    cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Y.real),
                    rstride=2, cstride=2)

    # Draw a set of x, y, z axes for reference.
    ax_lim = 0.5
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.
    ax.set_title(r'$Y_{{{},{}}}$'.format(el, m))
    ax_lim = 0.5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')

# plots all the SH terms in one figure
'''el_max = 3
figsize_px, DPI = 800, 100
figsize_in = figsize_px / DPI
fig = plt.figure(figsize=(figsize_in, figsize_in), dpi=DPI)
spec = gridspec.GridSpec(ncols=2*el_max+1, nrows=el_max+1, figure=fig)
for el in range(el_max+1):
    for m_el in range(-el, el+1):
        print(el, m_el)
        ax = fig.add_subplot(spec[el, m_el+el_max], projection='3d')
        #plot_Y(ax, el, m_el)  # plots in deta vs. dphi space
        #plot_Y_xyz(ax, el, m_el) # plots in cartesian space
plt.tight_layout()
plt.savefig('sph_harm.png')
plt.show()
exit()'''

fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(projection='3d')
# plot the sherical harmonics in the cartesian space
#plot_Y_xyz(ax, 1, 1)
#plt.savefig('Y{}_{}.png'.format(l, m))

fin = ROOT.TFile('deta_dphi.root')
h = fin.Get('h')
scale = h.Integral()
h.Scale(1/scale)
# plot spherical harmonics in delta eta vs. delta phi
#plot_Y(ax, 1, -1)
#plot data from file in the delta eta vs. delta phi space
plot_histo(ax, h)
#plot data from file in the cartesian space
#plot_histo_xyz(ax, h)

# code performing the SH parameters fit to the data
sh.SetParameter(0,1.0)
for ipar in range (1,16): sh.SetParameter(ipar, 0.0)
h.Fit(sh,'rn0')
#plot the sum of the SH terms fit to the data
#plot_sh(ax, sh)

# fill up a histogram with the SH parameters fit to the data
hpar = ROOT.TH1F('hpar','', 16, -0.5, 15.5)
for ipar in range(0, 16):
    hpar.SetBinContent(ipar+1, sh.GetParameter(ipar))
    hpar.SetBinError(ipar+1, sh.GetParError(ipar))

# plot the SH spectrum obtained from the fit
#figpars = plt.figure(figsize=plt.figaspect(1.))
#axpars = fig.add_subplot()
#plot_th1(hpar, axpars, color='black', markersize=10, markerfacecolor='gray')

plt.show()

