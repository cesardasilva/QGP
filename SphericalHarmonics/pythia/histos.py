import ROOT
import numpy as np
import matplotlib.pyplot as plt

class Hist(ROOT.TH1F):
    def __init__(self, par, name=''):
        super().__init__()
        if isinstance(par, np.ndarray):
            self.xbins = par
            super().SetName(name)
            super().SetBins(len(par)-1, np.asarray(par,'d'))
            super().Sumw2()
        elif isinstance(par, ROOT.TH1):
            super().SetName(par.GetName())
            bins = []
            nbins = par.GetNbinsX()
            for ibin in range(nbins):
                bins.append(par.GetXaxis().GetBinLowEdge(ibin+1))
            bins.append(par.GetXaxis().GetBinUpEdge(nbins))
            self.xbins = np.array(bins)
            super().SetBins(len(self.xbins)-1, np.asarray(self.xbins,'d'))
            super().Sumw2()        
            super().SetName(name)
            for ibin in range(nbins):
                super().SetBinContent(ibin+1, par.GetBinContent(ibin+1))
                super().SetBinError(ibin+1, par.GetBinError(ibin+1))
            super().SetEntries(par.GetEntries())
        else:
            print ('cannot identify initialization')
        
    def fill_array(self, xdata, w=None):
        nphisto, self.xbins = np.histogram(xdata, bins=self.xbins, weights=w)
        nphisto_err = np.sqrt(nphisto)
        for i in range(len(self.xbins)-1):
            val = super().GetBinContent(i+1) + nphisto[i]
            err = np.sqrt(super().GetBinContent(i+1) + nphisto_err[i])
            super().SetBinContent(i+1, val)
            super().SetBinError(i+1, err)
        super().SetEntries(len(xdata))

    def hist2array(self):
        npdata = np.zeros(len(self.xbins)-1)
        for i in range(len(self.xbins)-1):
            npdata[i] = super().GetBinContent(i+1)
        return npdata

    def get_binCenters(self):
        return 0.5*(xbins[1:] + xbins[:-1])

class Hist2D(ROOT.TH2F):
    def __init__(self, par, _ybins=None, name=None):
        if isinstance(par, np.ndarray):
            self.xbins = par
            self.ybins = _ybins
            super().__init__()
            super().SetName(name)
            super().SetBins(len(par)-1, np.asarray(par,'d'), len(_ybins)-1, np.asarray(_ybins,'d'))
            super().Sumw2()
        elif isinstance(par, ROOT.TH2): 
            super().__init__()
            binsX = []
            binsY = []
            nbinsX = par.GetNbinsX()
            nbinsY = par.GetNbinsY()
            for ibinY in range(nbinsY):
                binsY.append(par.GetYaxis().GetBinLowEdge(ibinY+1))
            binsY.append(par.GetYaxis().GetBinUpEdge(nbinsY))
            for ibinx in range(nbinsX):
                binsX.append(par.GetXaxis().GetBinLowEdge(ibinx+1))
            binsX.append(par.GetXaxis().GetBinUpEdge(nbinsX))
            self.xbins = np.array(binsX)
            self.ybins = np.array(binsY)
            super().SetBins(len(self.xbins)-1, self.xbins, len(self.ybins)-1, self.ybins)
            super().Sumw2()
            super().SetName(name)
            for ibiny in range(nbinsY):
                for ibinx in range(nbinsX):
                    super().SetBinContent(ibinx+1, ibiny+1, par.GetBinContent(ibinx+1, ibiny+1))
                    super().SetBinError(ibinx+1, ibiny+1, par.GetBinError(ibinx+1,ibiny+1))
            super().SetEntries(par.GetEntries())
        else:
            print ('cannot identify initialization')

           
    def fill_array(self, xdata, ydata, w=None):
        nphisto, self.xbins, self.ybins = np.histogram2d(xdata, ydata, bins=(self.xbins, self.ybins), weights=w)
        nphisto_noweights, self.xbins, self.ybins = np.histogram2d(xdata, ydata, bins=(self.xbins, self.ybins))
        rel_err = np.divide(1,nphisto_noweights, out=np.zeros_like(nphisto_noweights), where=nphisto_noweights!=0)
        nphisto_err = nphisto*np.sqrt(rel_err)
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                val = super().GetBinContent(i+1, j+1) + nphisto[i,j]
                err = np.sqrt(super().GetBinContent(i+1, j+1) + nphisto_err[i,j])
                super().SetBinContent(i+1, j+1, val)
                super().SetBinError(i+1, j+1, err)
        super().SetEntries(len(xdata))

    def copy_from(self, histo):
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                super().SetBinContent(i+1, j+1, histo.GetBinContent(i+1,j+1))
                super().SetBinError(i+1, j+1, histo.GetBinError(i+1,j+1))
        super().SetEntries(histo.GetEntries())      
        
    def hist2array(self):
        npdata = np.zeros((len(self.xbins), len(self.ybins)))
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                npdata[i,j] = super().GetBinContent(i+1, j+1)
        return npdata

    def get_binCenters(self):
        return (0.5*(xbins[1:] + xbins[:-1])), (0.5*(ybins[1:] + ybins[:-1])) 

    def matplotlib_draw(self,axs=None, vmin=None, vmax=None, cmap='rainbow'):
        hdata = self.hist2array()
        hdata = hdata.T
        min_cont = super().GetMinimum()
        if vmin!=None : min_cont = vmin
        max_cont = super().GetMaximum()
        if vmax!=None : max_cont = vmax
        X, Y = np.meshgrid(self.xbins, self.ybins)
        hdata[hdata==0] = None
        return axs.pcolormesh(X, Y, hdata, cmap=cmap, shading='auto',
                              vmin=min_cont, vmax=max_cont)
    

class Hist3D(ROOT.TH3F):
    def __init__(self, par, _ybins=None, _zbins=None, name=None):
        if isinstance(par, np.ndarray):
            self.xbins = par
            self.ybins = _ybins
            self.zbins = _zbins
            super().__init__()
            super().SetName(name)
            super().SetBins(len(par)-1, np.asarray(par,'d'),
                            len(_ybins)-1, np.asarray(_ybins,'d'),
                            len(_zbins)-1, np.asarray(_zbins,'d'))
            super().Sumw2()
        elif isinstance(par, ROOT.TH3): 
            super().__init__()
            binsX = []
            binsY = []
            binsZ = []
            nbinsX = par.GetNbinsX()
            nbinsY = par.GetNbinsY()
            nbinsZ = par.GetNbinsZ()
            for ibinZ in range(nbinsZ):
                binsZ.append(par.GetZaxis().GetBinLowEdge(ibinZ+1))
            binsZ.append(par.GetZaxis().GetBinUpEdge(nbinsZ))
            for ibinY in range(nbinsY):
                binsY.append(par.GetYaxis().GetBinLowEdge(ibinY+1))
            binsY.append(par.GetYaxis().GetBinUpEdge(nbinsY))
            for ibinX in range(nbinsX):
                binsX.append(par.GetXaxis().GetBinLowEdge(ibinX+1))
            binsX.append(par.GetXaxis().GetBinUpEdge(nbinsX))
            self.xbins = np.array(binsX)
            self.ybins = np.array(binsY)
            self.zbins = np.array(binsZ)
            super().SetName(name)
            super().SetBins(len(self.xbins)-1, self.xbins,
                            len(self.ybins)-1, self.ybins,
                            len(self.zbins)-1, self.zbins)
            super().Sumw2()
            for ibinz in range(nbinsZ):
                for ibiny in range(nbinsY):
                    for ibinx in range(nbinsX):
                        super().SetBinContent(ibinx+1, ibiny+1, ibinz+1,
                                              par.GetBinContent(ibinx+1, ibiny+1, ibinz+1))
                        super().SetBinError(ibinx+1, ibiny+1, ibinz+1,
                                            par.GetBinError(ibinx+1, ibiny+1, ibinz+1))
            super().SetEntries(par.GetEntries())
        else:
            print ('cannot identify initialization')
           
    def fill_array(self, xdata, ydata, zdata, w=None):
        nphisto, [self.xbins, self.ybins, self.zbins] = np.histogramdd((xdata, ydata, zdata), bins=(self.xbins, self.ybins, self.zbins), weights=w)
        nphisto_noweights, [self.xbins, self.ybins, self.zbins] = np.histogramdd((xdata, ydata, zdata), bins=(self.xbins, self.ybins, self.zbins))
        rel_err = np.divide(1,nphisto_noweights, out=np.zeros_like(nphisto_noweights), where=nphisto_noweights!=0)
        nphisto_err = nphisto*np.sqrt(rel_err)
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                for k in range(len(self.zbins)-1):
                    val = super().GetBinContent(i+1, j+1, k+1) + nphisto[i,j,k]
                    err = np.sqrt(super().GetBinContent(i+1, j+1, k+1) + nphisto_err[i,j,k])
                    super().SetBinContent(i+1, j+1, k+1, val)
                    super().SetBinError(i+1, j+1, k+1, err)
        super().SetEntries(len(xdata))

    def copy_from(self, histo):
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                for k in range(len(self.zbins)-1):
                    super().SetBinContent(i+1, j+1, k+1, histo.GetBinContent(i+1,j+1, k+1))
                    super().SetBinError(i+1, j+1, k+1, histo.GetBinError(i+1,j+1,k+1))
        super().SetEntries(histo.GetEntries())      
        
    def hist2array(self):
        npdata = np.zeros((len(self.xbins), len(self.ybins)))
        for i in range(len(self.xbins)-1):
            for j in range(len(self.ybins)-1):
                for k in range(len(self.zbins)-1):
                    npdata[i,j,k] = super().GetBinContent(i+1, j+1, k+1)
        return npdata

    def get_binCenters(self):
        return (0.5*(xbins[1:] + xbins[:-1])), (0.5*(ybins[1:] + ybins[:-1])), (0.5*(zbins[1:] + zbins[:-1]))


def plot_th1(histo, ax=None, xshift=0.0, xscale=1.0, show_negative=True, upperlim=False, is_box=False, **kargs):
    ax = ax if ax is not None else plt.gca()
    if 'markersize' not in kargs: kargs['markersize'] = 3
    if 'marker' not in kargs: kargs['marker'] = 'o'
    if 'ls' not in kargs: kargs['ls'] = 'none'
    if 'markeredgecolor' not in kargs: kargs['markeredgecolor'] = 'black'    
    x = []
    y = []
    xerr = []
    yerr = []
    xlim = []
    ylim = []
    xlimerr = []
    ylimerr = []    
    for i in range(1,histo.GetNbinsX()+1):
        if histo.GetBinError(i)==0: continue
        if np.isnan(histo.GetBinError(i)) : continue
        if (histo.GetBinContent(i)-histo.GetBinError(i)<0) and (show_negative==False):
            if upperlim :
                xlim.append((histo.GetXaxis().GetBinCenter(i)+xshift)*xscale)
                xlimerr.append(histo.GetXaxis().GetBinWidth(i)/2)
                ylim.append(histo.GetBinContent(i)+2*histo.GetBinError(i)) ## 90% CL
                ylimerr.append(histo.GetBinContent(i))
            continue
        x.append((histo.GetXaxis().GetBinCenter(i)+xshift)*xscale)
        xerr.append(histo.GetXaxis().GetBinWidth(i)/2)
        y.append(histo.GetBinContent(i))
        yerr.append(histo.GetBinError(i))
    if is_box == False:
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kargs)
    else :
        ax.bar(x, 2*np.array(yerr), bottom=np.array(y)-np.array(yerr), width=np.array(xerr), color=kargs['markerfacecolor'], edgecolor='black', align='center', alpha=kargs['alpha'])
    if len(xlim)>0 :
        ax.errorbar(xlim, ylim, xerr=xlimerr, yerr=ylimerr, color=kargs['markerfacecolor'],
                    uplims=True, ls='none')
    
def plot_th1_bars(histo, color=None, ax=None, **kargs):
    ax = ax if ax is not None else plt.gca()
    x = []
    ymin = []
    ymax = []
    for i in range(1,histo.GetNbinsX()+2):
        x.append(histo.GetXaxis().GetBinLowEdge(i))
        ymin.append(0)
        ymax.append(histo.GetBinContent(i)+histo.GetBinError(i))
    ax.fill_between(x, ymax, ymin, step='post', **kargs)
    
def errorfill(histo, color=None, ax=None, **kargs):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    x = []
    ymin = []
    ymax = []
    for i in range(1,histo.GetNbinsX()+1):
        x.append(histo.GetXaxis().GetBinCenter(i))
        ymin.append(histo.GetBinContent(i)-histo.GetBinError(i))
        ymax.append(histo.GetBinContent(i)+histo.GetBinError(i))
    ax.fill_between(x, ymax, ymin, color=color, **kargs) #label=histo.GetTitle())

def plot_tf1(func, ax=None, **kargs):
    ax = ax if ax is not None else plt.gca()
    x = []
    y = []
    step = (func.GetXmax()-func.GetXmin())/5000
    for i in range(0, 5000):
        x.append(func.GetXmin()+step*i)
        y.append(func.Eval(x[-1]))
    ax.plot(x, y, **kargs)

import csv

## https://tex.stackexchange.com/questions/180648/building-a-2d-histogram-with-pgfplots
def th2f_to_csv(hist, csv_file):
    """Print TH2F bin data to CSV file."""
    xbins, ybins = hist.GetNbinsX(), hist.GetNbinsY()
    xaxis, yaxis = hist.GetXaxis(), hist.GetYaxis()
    with open(csv_file, 'w') as f:
        c = csv.writer(f, delimiter=' ', lineterminator='\n')
        for ybin in xrange(1, ybins+2):
            y_lowedge = yaxis.GetBinLowEdge(ybin)
            for xbin in xrange(1, xbins+2):
                x_lowedge = xaxis.GetBinLowEdge(xbin)
                weight = hist.GetBinContent(xbin, ybin)
                c.writerow((x_lowedge, y_lowedge, weight))

## use in PGFplots
#\usepackage{tikz}
#\usepackage{pgfplots}
#\pgfplotsset{compat=newest}

#\begin{tikzpicture}
#  \begin{axis}[
#    view={0}{90},
#    colorbar,
#  ]
#    \addplot3[
#      surf,
#      shader=flat corner,
#      mesh/cols=51,
#      mesh/ordering=rowwise,
#    ] file {matrix.csv};
#  \end{axis}
#\end{tikzpicture}

