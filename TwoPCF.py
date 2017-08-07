import numpy as np
from matplotlib import pyplot as plt
import math as m
import sklearn
from sklearn.neighbors import BallTree
from astroML.correlation import two_point
from astroML.correlation import bootstrap_two_point


try:
    from sklearn.neighbors import KDTree
    sklearn_has_two_point = True
except ImportError:
    import warnings
    sklearn_has_two_point = False


def two_point(data, bins, method='standard', data_R=None, random_state=None):

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")
    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")
    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    print 'ndim'
    print data.ndim
    n_samples, n_features = data.shape

    Nbins = len(bins) - 1

    # shuffle all but one axis to get background distribution

    factor = len(data_R) * 1. / len(data)
    if sklearn_has_two_point:
        # Fast two-point correlation functions added in scikit-learn v. 0.14
        KDT_D = KDTree(data)
        KDT_R = KDTree(data_R)
        counts_DD = KDT_D.two_point_correlation(data, bins)
        counts_RR = KDT_R.two_point_correlation(data_R, bins)

    else:
        warnings.warn("Version 0.3 of astroML will require scikit-learn "
                      "version 0.14 or higher for correlation function "
                      "calculations. Upgrade to sklearn 0.14+ now for much "
                      "faster correlation function calculations.")

        BT_D = BallTree(data)
        BT_R = BallTree(data_R)

        counts_DD = np.zeros(Nbins + 1)
        counts_RR = np.zeros(Nbins + 1)

        for i in range(Nbins + 1):
            counts_DD[i] = np.sum(BT_D.query_radius(data, bins[i], count_only=True))
            counts_RR[i] = np.sum(BT_R.query_radius(data_R, bins[i], count_only=True))

    DD1 = np.diff(counts_DD)
    RR1 = np.diff(counts_RR)

    RR_zero = (RR1 == 0)

    RR1[RR_zero] = 1

    if sklearn_has_two_point:
        counts_DR = KDT_R.two_point_correlation(data, bins)
    else:
        counts_DR = np.zeros(Nbins + 1)
        for i in range(Nbins + 1):
            counts_DR[i] = np.sum(BT_R.query_radius(data, bins[i], count_only=True))

    DR1 = np.diff(counts_DR)
    corr1 = (factor ** 2 * DD1 - 2 * factor * DR1 + RR1) / RR1
    corr1[RR_zero] = np.nan
    return  corr1


