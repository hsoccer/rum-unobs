import numpy as np
import pandas as pd

from .utils import get_power_sets


class LatticeElement():
    """
    class for subsets
    """
    def __init__(self, D):
        self.D = D


def get_lattice_instances(X):
    """
    Parameters
    ----------
    X : iterable
        universal set
    
    Returns
    -------
    list_insts : list
        list of instances of each subset
    """
    # X = (1, 2, 3)
    all_subsets = list(map(set, get_power_sets(X)))
    list_insts = []
    for D in all_subsets:
        # D = all_subsets[2]
        inst = LatticeElement(tuple(D))
        inst.supsets = [tuple(subset) for subset in all_subsets if (subset != D) & (D.issubset(subset))]
        inst.subsets = [tuple(subset) for subset in all_subsets if not (D.issubset(subset))]
        list_insts.append(inst)
    return list_insts


def BM(rho, X, add_trivial=True):
    """
    Compute all the BM polynomials of rho for each choice set and element in it.

    Parameters
    ----------
    rho : function
        a random choice function, which takes (D, x) as its arguments
    X : iterable
        universal set
    add_trivial : boolian (default = True)
        if true, add singleton choice sets

    Returns
    ------- 
    df_bm : pd.DataFrame
        index : D / columns : x

    Examples
    --------
    >>> X = (1, 2, 3)
    >>> rho_uniform = lambda D, x: 1 / len(D)
    >>> rho = rho_uniform
    >>> BM(rho, X)
    x                 1         2         3
    D                                      
    (1,)       0.333333       NaN       NaN
    (2,)            NaN  0.333333       NaN
    (3,)            NaN       NaN  0.333333
    (1, 2)     0.166667  0.166667       NaN
    (1, 3)     0.166667       NaN  0.166667
    (2, 3)          NaN  0.166667  0.166667
    (1, 2, 3)  0.333333  0.333333  0.333333

    Notes
    -----
    When some choice data are missing in rho, some BM polynomials are not computable.
    In this case, such cells are left blank.
    """
    if add_trivial:
        def rho_added(D, x):
            if (len(D) == 1) and (x in D):
                return 1
            else:
                return rho(D, x)
    else:
        rho_added = rho

    X = tuple(X)
    list_lattice_instances = get_lattice_instances(X)
    bm = dict()
    for inst in list_lattice_instances[::-1]:
        D = inst.D
        if D == X:
            for x in D:
                try:
                    bm[(D, x)] = rho_added(D, x)  # rho(D, x) is available
                except:
                    bm[(D, x)] = np.nan  # rho(D, x) is not available
        else:
            all_supsets = inst.supsets
            for x in D:
                try:
                    bm[(D, x)] = rho_added(D, x) - sum([bm[(supset, x)] for supset in all_supsets])  # rho(D, x) is available
                except:
                    bm[(D, x)] = np.nan  # rho(D, x) is not available
    df_bm_stacked = pd.DataFrame([(k[0], k[1], v) for k, v in bm.items()])
    df_bm_stacked.columns = ["D", "x", "value"]
    df_bm = df_bm_stacked.pivot(index="D", columns="x", values="value").sort_index(key=lambda lst: [(len(elem), elem) for elem in lst])
    return df_bm


def observable_net_outflow(rho, X, add_trivial=True):
    """
    Compute resource constraints on the network.
    
    Parameters
    ----------
    rho : function
        a random choice function, which takes (D, x) as its arguments
    X : iterable
        universal set

    Returns
    ------- 
    delta : dict
        index : D / value : net outflow    
    """
    bm = BM(rho, X, add_trivial)
    
    X = tuple(X)
    list_lattice_instances = get_lattice_instances(X)
    
    delta = dict()
    
    for inst in list_lattice_instances:
        D = inst.D
        if D == X:
            delta[D] = 1 - bm.loc[[D], :].sum().sum()
        elif D == ():
            delta[D] = sum([bm.loc[[(x,)], x].sum() for x in X]) - 1
        else:
            inflow = bm.loc[[D], :].sum().sum()
            outflow = 0
            for x in set(X).difference(D):
                E = tuple(sorted(set(D) | {x}))
                outflow += bm.loc[[E], x].sum()
            delta[D] = outflow - inflow
    
    return delta

