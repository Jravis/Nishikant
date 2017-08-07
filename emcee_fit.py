import numpy as np
import matplotlib.pyplot as plt
import lmfit
import emcee
import corner

name = '/dataspace/sandeep/Nishi_plots/for_sandeep_saili/modes_037.dat'

Delta_m_r = np.genfromtxt(name, usecols=0)
Delta_m_i = np.genfromtxt(name, usecols=1)
Delta_f_r = np.genfromtxt(name, usecols=2)
Delta_f_i = np.genfromtxt(name, usecols=3)
mu_sqr= np.genfromtxt(name, usecols=5)

plt.plot(Delta_m_r, Delta_f_r, 'b.')
plt.plot(Delta_m_r, -0.12*(1. + 1.47*mu_sqr)*Delta_m_r, 'r*', ms=6)

"""
Delta_m = []
Delta_f = []
Mu_sqr = []
for ii in xrange(len(mu_sqr)):
    Delta_m.append(Delta_m_r[ii])
    Delta_m.append(Delta_m_i[ii])
    Delta_f.append(Delta_f_r[ii])
    Delta_f.append(Delta_f_i[ii])
    Mu_sqr.append(mu_sqr[ii])
    Mu_sqr.append(mu_sqr[ii])

Delta_f = np.asarray(Delta_f, dtype=np.float64)
Delta_m = np.asarray(Delta_m, dtype=np.float64)
Mu_sqr= np.asarray(Mu_sqr, dtype=np.float64)
"""

pars = lmfit.Parameters()
"""
pars.add('beta', value=2.0, min=0.05, max=10, brute_step=0.05)
pars.add('bias', value=-2.0, min=-5.0, max=2.0, brute_step=0.03517588)
pars.add('lnPn', value = 5, min=-10,max= 10,  brute_step= 0.10050251)
"""


"""
pars.add('beta', value=2.0, min=0.05, max=10)
pars.add('bias', value=-2.0, min=-5.0, max=2.0)
pars.add('lnPn', value = 5, min=-10,max= 10)



def residual(p1):    
    
    beta = p1['beta'].value
    bias = p1['bias'].value
    lnPn = p1['lnPn'].value
    
    
    tmp = bias*(1+beta*mu_sqr)
    tmp = Delta_f_r-(Delta_m_r*tmp)

    tmp1 = bias*(1. + beta*mu_sqr)
    tmp1 = Delta_f_i-(Delta_m_i*tmp1)

    return  (tmp*tmp)/np.exp(lnPn) + np.log(2.*np.pi)*lnPn + (tmp1*tmp1)/np.exp(lnPn) + np.log(2.*np.pi)*lnPn 




def lnprob(p1):
    resid = residual(p1)
    return -0.5 * np.sum(resid)



mi = lmfit.Minimizer(residual, pars)
#out1 = mi.minimize(method='brute')
out1 = mi.minimize(method='differential_evolution')
#out1 = mi.minimize(method='Nelder')

#print out1.params

mini = lmfit.Minimizer(lnprob, out1.params)
#mini = lmfit.Minimizer(lnprob, pars)
res = mini.emcee(burn=300, steps=1000, thin=10, nwalkers= 100,  is_weighted=False, params=mi.params, seed=1230127)



print('------------------------------------------')
print("median of posterior probability distribution")
print('------------------------------------------')
lmfit.report_fit(res.params)

fig1 = corner.corner(res.flatchain, bins=100, labels=res.var_names,
        truths=list(res.params.valuesdict().values()), 
        quantiles=[0.16,0.5, 0.84], show_titles=True, labels_args={"fontsize": 40})
                    
fig1.set_size_inches(10,10)
fig1.savefig('/dataspace/sandeep/Nishi_plots/for_sandeep_saili/lymanalfa_forest_Scatter_10000step.pdf', dpi=300)

highest_prob = np.argmax(res.lnprob)
hp_loc = np.unravel_index(highest_prob, res.lnprob.shape)
mle_soln = res.chain[hp_loc]
for i, par in enumerate(pars):
    pars[par].value = mle_soln[i]

print("\nMaximum likelihood Estimation")
print('-----------------------------')
print(pars)



print("\nFinally lets work out 2-sigma error estimate")

quantiles = np.percentile(res.flatchain['beta'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in a", 0.5 * (quantiles[-1] - quantiles[0]))

quantiles = np.percentile(res.flatchain['bias'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in b", 0.5 * (quantiles[-1] - quantiles[0]))

quantiles = np.percentile(res.flatchain['lnPn'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in c", 0.5 * (quantiles[-1] - quantiles[0]))



plt.figure(2, figsize=(8, 6))
plt.plot(res.flatchain['lnPn'], 'k-', alpha=0.3)
#plt.savefig('/dataspace/sandeep/Nishi_plots/for_sandeep_saili/lnPn_chain.pdf', dpi=300)

plt.figure(3, figsize=(8, 6))
plt.plot(res.flatchain['bias'], 'r-', alpha=0.3)
#plt.savefig('/dataspace/sandeep/Nishi_plots/for_sandeep_saili/bias_chain.pdf', dpi=300)

plt.figure(4, figsize=(8, 6))
plt.plot(res.flatchain['beta'], 'b-', alpha=0.3)
#plt.savefig('/dataspace/sandeep/Nishi_plots/for_sandeep_saili/beta_chain.pdf', dpi=300)
"""



plt.show()
