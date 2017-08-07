import numpy as np
import matplotlib.pyplot as plt
import lmfit
import emcee
import corner

name = "/dataspace/sandeep/Nishi_plots/Test_emceedata/data_cubic_sample_0.dat"
x = np.genfromtxt(name, usecols = 0) 
y = np.genfromtxt(name, usecols = 1) 
y_err= np.genfromtxt(name, usecols = 2) 
data = np.zeros((3, len(x)), dtype=np.double)
data[0, :] = x
data[1, :] = y
data[2, :] = y_err


pars = lmfit.Parameters()
pars.add_many( ('a', 5*10**-4), ('b', 0.001), ('c',1.0 ), ('d',-0.001) ) 


def residual(p1):
    a = p1['a'].value
    b = p1['b'].value
    c = p1['c'].value
    d = p1['d'].value
    model = a*data[0, :]**3 + b*data[0,:]**2+c*data[0,:]+d
    return data[1, :] -  model 


def lnprob(p):
    resid = residual(p)
    s = data[2, :]
    resid *= 1 / s
    resid *= resid
    resid += np.log(2 * np.pi * s ** 2)
    return -0.5 * np.sum(resid)

mi = lmfit.Minimizer(residual, pars)
# first solve with Nelder-Mead
out1 = mi.minimize(method='Nelder')


print('------------------------------------------')

mini = lmfit.Minimizer(lnprob, out1.params)
res = mini.emcee(burn=400, steps=2000, thin=10)#, params=mi.params)

print("median of posterior probability distribution")
print('------------------------------------------')

lmfit.report_fit(res.params)



fig1 = corner.corner(res.flatchain, labels=res.var_names, truths=list(res.params.valuesdict().values()),
        truth_color='skyblue')

fig1.savefig('/dataspace/sandeep/Nishi_plots/Test_emceedata/emcee_sample0_fit.eps', dpi=100)


highest_prob = np.argmax(res.lnprob)
hp_loc = np.unravel_index(highest_prob, res.lnprob.shape)
mle_soln = res.chain[hp_loc]




for i, par in enumerate(pars):
    pars[par].value = mle_soln[i]

print("\nMaximum likelihood Estimation")
print('-----------------------------')
print(pars)


plt.figure(2, figsize=(8, 6))
plt.plot(res.flatchain['a'], 'k-', alpha=0.3)


# Finally lets work out a 1 and 2-sigma error estimate for 't1'

quantiles = np.percentile(res.flatchain['a'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in a", 0.5 * (quantiles[-1] - quantiles[0]))
a_2sig = 0.5 * (quantiles[-1] - quantiles[0])


quantiles = np.percentile(res.flatchain['b'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in b", 0.5 * (quantiles[-1] - quantiles[0]))
b_2sig = 0.5 * (quantiles[-1] - quantiles[0])


quantiles = np.percentile(res.flatchain['c'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in c", 0.5 * (quantiles[-1] - quantiles[0]))
c_2sig = 0.5 * (quantiles[-1] - quantiles[0])

quantiles = np.percentile(res.flatchain['d'], [2.28, 15.9, 50, 84.2, 97.7])
print("2 sigma spread in  d", 0.5 * (quantiles[-1] - quantiles[0]))
d_2sig = 0.5 * (quantiles[-1] - quantiles[0])

a = 0.973
b = -10.0755
c = -3.5187
d = 4.2621
model= a*x**3 + b*x**2+c*x+d

err = (0.045520)*x**3 + (0.177237)*x**2 + (0.301567)*x + (0.227867)

fig1 = plt.figure(3, figsize=(8, 6))
plt.errorbar(x, y, yerr=y_err, fmt='o', label='sample0')
plt.errorbar(x, model, yerr=err, fmt='o', label='emcee')
plt.legend()
#plt.fill_between(x, model-plus_model, model+plus_model,
#        facecolor='lightgrey', edgecolors="k")

fig1.savefig('sample0.pdf', dpi=300)




plt.show()
