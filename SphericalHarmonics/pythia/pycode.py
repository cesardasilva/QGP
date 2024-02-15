
import sys
import numpy as np
from math import exp, pi
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
import ROOT
plt.rc('text', usetex=False)
from histos import Hist2D, Hist3D, plot_th1, errorfill, plot_tf1

#l, m = int(sys.argv[1]), int(sys.argv[2])

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
    theta = x[1]
    total = 0.0
    ipar = 0
    for l in range(0,4):
        for m in range(-l,l+1):
            Y = sph_harm(abs(m), l, phi, theta)
            if m < 0: Y = np.sqrt(2) * (-1)**m * Y.imag
            elif m > 0: Y = np.sqrt(2) * (-1)**m * Y.real
            total += pars[ipar]*Y
            ipar += 1
    return total

sh = ROOT.TF2('sh', SH, -pi, pi, -3, 3, 16)

# Grids of polar and azimuthal angles
thetabins = np.linspace(-3, 3, 50)
phi = np.linspace(-pi, pi, 50)
# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(thetabins, phi)
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(theta) * np.sin(phi),
                np.sin(theta) * np.cos(phi),
                np.cos(theta)])

def plot_histo(ax,h) :
    scale = h.Integral()
    h.Scale(1/scale)
    data = []
    for x in range(1,h.GetNbinsX()+1):
        datay = [] #np.zeros(h.GetNbinsY()+1, dtype='float32')
        for y in range(1,h.GetNbinsY()+1):
            #eta = h.GetBinCenter(y)
            #theta = abs(2*exp(-eta))
            #thetabin = np.digitize(theta, thetabins) 
            datay.append(h.GetBinContent(x,y))
        data.append(datay)
    data = np.array(data)
    Yx, Yy, Yz = np.abs(data) * xyz
    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('PRGn'))
    cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Yz),
                    rstride=2, cstride=2)    
    
def plot_Y(ax, el, m):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""

    # NB In SciPy's sph_harm function the azimuthal coordinate, theta,
    # comes before the polar coordinate, phi.
    '''Y = sph_harm(abs(m), el, phi, theta)

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    Yx, Yy, Yz = np.abs(Y) * xyz'''

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

'''el_max = 3
figsize_px, DPI = 800, 100
figsize_in = figsize_px / DPI
fig = plt.figure(figsize=(figsize_in, figsize_in), dpi=DPI)
spec = gridspec.GridSpec(ncols=2*el_max+1, nrows=el_max+1, figure=fig)
for el in range(el_max+1):
    for m_el in range(-el, el+1):
        print(el, m_el)
        ax = fig.add_subplot(spec[el, m_el+el_max], projection='3d')
        plot_Y(ax, el, m_el)
plt.tight_layout()
plt.savefig('sph_harm.png')
plt.show()'''

fin = ROOT.TFile('deta_dphi.root')
h = fin.Get('h')
scale = h.Integral()
h.Scale(1/scale)
sh.SetParameter(0,1.0)
for ipar in range (1,16): sh.SetParameter(ipar, 0.0)
h.Fit(sh,'rn0')

hpar = ROOT.TH1F('hpar','', 16, -0.5, 15.5)
for ipar in range(0, 16):
    hpar.SetBinContent(ipar+1, sh.GetParameter(ipar))
    hpar.SetBinError(ipar+1, sh.GetParError(ipar))

fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(projection='3d')
plot_histo(ax, h)
plot_Y(ax, l, m)
plt.savefig('Y{}_{}.png'.format(l, m))

figpars = plt.figure(figsize=plt.figaspect(1.))
axpars = fig.add_subplot()

plot_th1(hpar, axpars, color='black', markersize=10, markerfacecolor='gray')

plt.show()

