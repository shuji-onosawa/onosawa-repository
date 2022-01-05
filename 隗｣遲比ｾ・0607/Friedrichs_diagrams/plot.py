import sympy as sp
from sympy.plotting import plot
import numpy as np

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# setup some plot defaults
linewidth = 2
fontsize  = 30
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Times')
plt.rc('font', size=fontsize)
plt.rc('axes', linewidth=linewidth)
plt.rc('axes', labelsize=fontsize)
plt.rc('legend', fontsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('xtick', top=True)
plt.rc('xtick.major', width=linewidth)
plt.rc('xtick.major', size=15)
plt.rc('xtick.minor', width=linewidth)
plt.rc('xtick.minor', size=8)
plt.rc('ytick', labelsize=fontsize)
plt.rc('ytick', right=True)
plt.rc('ytick.major', width=linewidth)
plt.rc('ytick.major', size=15)
plt.rc('ytick.minor', width=linewidth)
plt.rc('ytick.minor', size=8)
# plt.rc('xtick', labelbottom='off')
plt.rc('xtick', direction='in')
# plt.rc('ytick', labelleft='off')
plt.rc('ytick', direction='in')
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble  


# Calculate non-relativistic group velocity >>>>>>>>>>>>>>>>>>
def vgr_NR(vph, Cs, Va, t):
  cost = np.cos(t)
  n = np.array([np.sin(t), np.cos(t)])
  b = np.array([0.0, 1.0])
  denominator = 3.0*vph**5 - 2.0*( Cs**2 + (1.0 + cost**2)*Va**2 )*vph**3 + Va**2*(2.0*Cs**2 + Va**2)*cost**2*vph \
               - kdi**2*Va**2*cost**2*(2.0*vph**2 - Cs**2)*vph
  numerator   =  ( (Cs**2 + Va**2)*vph**4 - Va**2*(2.0*Cs**2 + Va**2)*cost**2*vph**2 + Cs**2*Va**4*cost**4 )*n \
               + ( Va**2*cost*vph**4 - Va**2*(2.0*Cs**2 + Va**2)*cost*vph**2 + 2.0*Cs**2*Va**4*cost**3  )*b \
               + kdi**2*Va**2*(\
                    (-2.0*Cs**2 + vph**2)*cost**2*vph**2*n \
                  + (vph**2 - Cs**2)*cost*vph**2*b \
                 )

  vgr_perp, vgr_para = numerator/denominator

  return [vgr_para, vgr_perp]
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

betas = [0.5, 1.0, 5.0]
gamma = 5.0/3.0

for beta in betas:
  alp   = beta*gamma/2.0 # = cS^2 / vA^2
  sgm   = 4.0*alp/(1.0 + alp)**2

  epsilon = 1e-3
  theta = np.linspace(0.0 + epsilon, 2*np.pi - epsilon, 100)

  def n_vec(t):
    return np.asarray([np.cos(t), np.sin(t)])

  def t_vec(t):
    return np.asarray([-np.sin(t), np.cos(t)])

  def vph_a(t):
    return abs(np.cos(t))

  def vph_s(t):
    return np.sqrt(0.5*(1.0 + alp - np.sqrt((1.0 + alp)**2 - 4.0*alp*(np.cos(t))**2)))

  def vph_f(t):
    return np.sqrt(0.5*(1.0 + alp + np.sqrt((1.0 + alp)**2 - 4.0*alp*(np.cos(t))**2)))

  def vgr_s_corr(t):
    cost = np.cos(t)
    sint = np.sin(t)
    return sgm*sint*cost/(2.0*np.sqrt(1.0 - sgm*cost**2)*(1 - np.sqrt(1.0 - sgm*cost**2)))

  def vgr_f_corr(t):
    cost = np.cos(t)
    sint = np.sin(t)
    return sgm*sint*cost/(2.0*np.sqrt(1.0 - sgm*cost**2)*(1 + np.sqrt(1.0 - sgm*cost**2)))

  vph_a_vec = np.asarray([vph_a(t)*n_vec(t) for t in theta])
  vph_s_vec = np.asarray([vph_s(t)*n_vec(t) for t in theta])
  vph_f_vec = np.asarray([vph_f(t)*n_vec(t) for t in theta])

  vgr_a_vec = np.asarray([[np.sign(np.sin(t)), 0.0] for t in theta])
  vgr_s_vec = np.asarray([vph_s(t)*(n_vec(t) - vgr_s_corr(t)*t_vec(t)) for t in theta])
  vgr_f_vec = np.asarray([vph_f(t)*(n_vec(t) + vgr_f_corr(t)*t_vec(t)) for t in theta])

  fig, ax = plt.subplots()
  plt.figure(figsize=(14,8))

  plt.subplot(121,aspect=1.0)
  plt.plot( vph_f_vec.T[0], vph_f_vec.T[1], color=plt.get_cmap('tab10').colors[2], lw=2, label=r'Fast')
  plt.plot( vph_s_vec.T[0], vph_s_vec.T[1], color=plt.get_cmap('tab10').colors[1], lw=2, label=r'Slow')
  plt.plot( vph_a_vec.T[0], vph_a_vec.T[1], color=plt.get_cmap('tab10').colors[0], lw=2, label=r'Alfv\'{e}n')
  plt.xlabel(r'$v^{\mr{ph}}_z/v_\rmA$')
  plt.ylabel(r'$v^{\mr{ph}}_x/v_\rmA$')
  legend = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=25)
  legend.set_title(r'$\beta$ = %.1f' % beta, prop={'size': 25})

  plt.subplot(122,aspect=1.0)
  plt.plot( vgr_f_vec.T[0], vgr_f_vec.T[1]     , color=plt.get_cmap('tab10').colors[2], lw=2, label=r'Fast')
  plt.plot( vgr_s_vec.T[0], vgr_s_vec.T[1]     , color=plt.get_cmap('tab10').colors[1], lw=2, label=r'Slow')
  plt.plot( vgr_a_vec.T[0], vgr_a_vec.T[1], 'o', color=plt.get_cmap('tab10').colors[0], ms=5, label=r'Alfv\'{e}n')
  plt.xlabel(r'$v^{\mr{gr}}_z/v_\rmA$')
  plt.ylabel(r'$v^{\mr{gr}}_x/v_\rmA$')

  plt.tight_layout()
  plt.savefig('beta%.2f.pdf' % beta)
