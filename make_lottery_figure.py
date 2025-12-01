import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.optimize import minimize
from itertools import permutations, combinations

from packages.utils import from_df_to_func, get_power_sets
from packages.bm import BM, observable_net_outflow

plt.style.use("ggplot")
plt.rcParams["font.family"] = "DejaVu Serif"


# fitting data to logit
def get_choice_prob(X, D_cal, list_utility):
    """
    Return the choice probability vector implied by the logit model based on LIST_UTILITY.
    """
    dict_utility_exp = dict(zip(X, [np.exp(utility) for utility in list_utility]))
    list_choice_prob = []
    for D in D_cal:
        denom = sum([dict_utility_exp[x] for x in D])
        choice_prob = [dict_utility_exp[x] / denom if x in D else 0 for x in X]
        list_choice_prob.append(choice_prob)
    df_choice_prob = pd.DataFrame(list_choice_prob, columns = X, index = D_cal)
    return df_choice_prob


def get_resid(df_choice_prob_true, list_utility_wo_0):
    X = tuple(df_choice_prob_true.columns)
    D_cal = list(df_choice_prob_true.index)
    list_utility = [0] + list(list_utility_wo_0)
    df_choice_prob = get_choice_prob(X, D_cal, list_utility)
    resid = np.linalg.norm((df_choice_prob_true - df_choice_prob).stack())
    return resid


def get_calibrated_data(rho_df):
    """
    Find the rationalizable dataset closest to RHO_DF.
    """
    x0 = np.random.uniform(low = -3, high = 3, size = 4)
    res = minimize(lambda x: get_resid(rho_df, x), x0, method="nelder-mead")
    list_utility = [0] + list(res.x)
    rho_fit_df = get_choice_prob(X, list(rho_df.index), list_utility)
    return rho_fit_df


def from_mu_to_choice_prob(mu, prefs):
    """
    Args:
        mu (np.array): probability measures over PREFS
        prefs (list)
    """
    mu = np.hstack([mu, [1 - mu.sum()]])
    X = tuple(sorted(prefs[0]))
    rows = []
    for i in range(2, len(X) + 1):
        rows += list(combinations(X, i))
    cols = X
    rho = np.zeros((len(rows), len(cols)))
    
    for row, D in enumerate(rows):
        for col, x in enumerate(X):
            # row, D = 0, (0, 1)
            # col, x = 0, 0
            if x not in D:
                rho[row, col] = np.nan
            else:
                rho[row, col] = np.sum((np.array([[y for y in pref if y in D][0] for pref in prefs]) == x) * mu)
                
    return pd.DataFrame(rho, columns=cols, index=rows) 


def get_calibrated_data_nonparametric(rho_df):
    X = rho_df.columns.tolist()
    prefs = list(permutations(X))  # (0, 1, 2, 3, 4) means 0 > 1 > 2 > 3 > 4
    num_prefs = len(prefs)
    
    mu0 = np.ones(num_prefs - 1) / num_prefs
    res = minimize(
            lambda mu: ((from_mu_to_choice_prob(mu, prefs) - rho_df)**2).sum().sum(), 
            mu0,
            bounds=[(0, None) for _ in range(num_prefs - 1)],
            constraints=({'type': 'ineq', 'fun': lambda mu: 1 - mu.sum()},),
        )
    rho_hat_df = from_mu_to_choice_prob(res.x, prefs)
    return rho_hat_df


# 0 & 1 are not observed + cut-off
def get_masked_data(rho_df, X_unobs):
    X = tuple(rho_df.columns)
    rho_missing_df = rho_df.copy()
    rho_missing_df.iloc[:, X_unobs] = np.nan
    
    # D contains no unobservables
    for idx, D in enumerate(rho_missing_df.index):
        for x in set(X).difference(D):
            rho_missing_df.iloc[idx, x] = 0
            
    # D contains only one unobservable
    for idx, D in enumerate(rho_missing_df.index):
        if rho_missing_df.iloc[idx, :].isna().sum() == 1:
            rho_missing_df.iloc[idx, :] = rho_missing_df.iloc[idx, :].fillna(1 - rho_missing_df.iloc[idx, :].sum())

    return rho_missing_df


def get_prob_bounds(rho_fitted_df, X_obs, X_unobs, x_obj):
    
    # not observed + cut-off
    rho_missing_df = get_masked_data(rho_fitted_df, X_unobs)
    rho_missing_func = from_df_to_func(rho_missing_df)
    bm_missing = BM(rho_missing_func, X)
    
    # compute delta for each choice set
    net_outflow_dict = observable_net_outflow(rho_missing_func, X)

    # all edges
    edge_list = []
    for D in get_power_sets(X):
        for x in D:
            E = tuple(sorted(set(D).difference({x})))
            edge_list.append((E, D))

    # (un)observable edges
    # change based on missing structure    
    starts = list(bm_missing.stack().reset_index().apply(lambda df: tuple(sorted(set(df["D"]).difference({df["x"]}))), axis = 1))
    ends = list(bm_missing.stack().reset_index()["D"])
    obs_edge_list = [(E, D) for E, D in zip(starts, ends)]
    unobs_edge_list = [(E, D) for E, D in edge_list if not (E, D) in obs_edge_list]

    # make incidence matrix for unobservable edges
    incidence_matrix_list = []
    for edge in unobs_edge_list:
        E, D = edge
        incidence_matrix_list.append([edge, E, -1])
        incidence_matrix_list.append([edge, D, 1])
    incidence_matrix_stack = pd.DataFrame(incidence_matrix_list, columns=["edge", "vertex", "value"])
    incidence_matrix = incidence_matrix_stack.pivot(index="vertex", columns="edge", values="value").sort_index(key=lambda lst: [(len(elem), elem) for elem in lst]).astype("Int64").fillna(0)
    incidence_matrix = incidence_matrix[pd.Index(unobs_edge_list)]

    # linear programming for bounds
    # compare identified sets
    A = incidence_matrix.values
    b = pd.DataFrame(net_outflow_dict.items()).set_index(0).values

    num_vertex, num_edge = A.shape
    
    target_choice_sets = [D for D in rho_missing_df[[all([x in idx for x in X_unobs])for idx in rho_missing_df.index]].index if x_obj in D]
    target_indices = [idx for idx, D in enumerate(rho_missing_df.index) if D in target_choice_sets]
    naive_upper_df = 1 - rho_missing_df.iloc[target_indices, :].loc[:, X_obs].sum(axis=1)  # upper bounds should be based on observed rather than on fitted
    L_sum = rho_missing_df.iloc[ :].loc[:, X_obs].sum(axis=1) 

    id_set_list = []
    for D in naive_upper_df.index:
        naive_upper = naive_upper_df[D] 
        c_lower = np.zeros(num_edge)
        for E in get_power_sets(X):
            if set(D).issubset(E):
                edge = (tuple(sorted(set(E).difference({x_obj}))), E)
                idx = unobs_edge_list.index(edge)
                c_lower[idx] = 1
        c_upper = - c_lower
        rum_lower = linprog(c_lower, A_ub=A, b_ub=b).fun
        rum_upper = -linprog(c_upper, A_ub=A, b_ub=b).fun
        if x_obj == 0:
            if D != (0,1):
                mon_lower = L_sum[D[1:]] - L_sum[D]
                mon_upper = 1- L_sum[D[:1] + D[2:]]
            else:
                mon_lower = 0
                mon_upper = 1
        if x_obj == 1:
            if D != (0,1):
                mon_lower = L_sum[D[:1] + D[2:]] - L_sum[D]
                mon_upper = 1- L_sum[D[1:]]
            else:
                mon_lower = 0
                mon_upper = 1
        id_set_list.append([D, rum_lower, rum_upper, 0, naive_upper, mon_lower, mon_upper])
    id_set_df = pd.DataFrame(id_set_list, columns=["D", "RUM LB", "RUM UB", "Naive LB", "Naive UB", "Mon LB", "Mon UB"]).set_index("D")
    
    return id_set_df



if __name__ == "__main__":
    rho_df = pd.read_csv("collective_choice_data.csv", index_col=0)
    rho_df.index = [eval(elem) for elem in rho_df.index]
    rho_df = rho_df.sort_index(key=lambda lst: [(len(elem), elem) for elem in lst])
    rho_df.columns = [eval(elem) for elem in rho_df.columns]

    X = tuple(rho_df.columns)
    X_unobs = (0, 1)
    x_obj = 0
    X_obs = tuple(sorted(set(X).difference(X_unobs)))
    
    # calibrate data
    # rho_fitted_df = get_calibrated_data(rho_df)
    rho_fitted_df = get_calibrated_data_nonparametric(rho_df)
    
    # get bounds
    id_set_df_0 = get_prob_bounds(rho_fitted_df, X_obs, X_unobs, x_obj=0)
    id_set_df_1 = get_prob_bounds(rho_fitted_df, X_obs, X_unobs, x_obj=1)
    target_choice_sets = id_set_df_0.index
    
    # visualize
    X_axis = np.arange(len(target_choice_sets))

    # error bar
    rum_bdd_mid_array_0 = (id_set_df_0.loc[:, "RUM UB"] + id_set_df_0.loc[:, "RUM LB"]).values / 2
    rum_error_array_0 = (id_set_df_0.loc[:, "RUM UB"] - id_set_df_0.loc[:, "RUM LB"]).values / 2
    rum_bdd_mid_array_1 = (id_set_df_1.loc[:, "RUM UB"] + id_set_df_1.loc[:, "RUM LB"]).values / 2
    rum_error_array_1 = (id_set_df_1.loc[:, "RUM UB"] - id_set_df_1.loc[:, "RUM LB"]).values / 2
    
    mon_bdd_mid_array_1 = (id_set_df_1.loc[:, "Mon UB"] + id_set_df_1.loc[:, "Mon LB"]).values / 2
    mon_error_array_1 = (id_set_df_1.loc[:, "Mon UB"] - id_set_df_1.loc[:, "Mon LB"]).values / 2
    mon_bdd_mid_array_0 = (id_set_df_0.loc[:, "Mon UB"] + id_set_df_0.loc[:, "Mon LB"]).values / 2
    mon_error_array_0 = (id_set_df_0.loc[:, "Mon UB"] - id_set_df_0.loc[:, "Mon LB"]).values / 2

    naive_bdd_mid_array = (id_set_df_0.loc[:, "Naive UB"] + id_set_df_0.loc[:, "Naive LB"]).values / 2
    naive_error_array = (id_set_df_0.loc[:, "Naive UB"] - id_set_df_0.loc[:, "Naive LB"]).values / 2

    # plot
    plt.figure(figsize=(10, 6))

    plt.scatter(X_axis - 0.2, rho_fitted_df.loc[target_choice_sets, 0], marker='o', color="red", facecolors='none', label=r"Fitted value ($x = 0$)", s=50)
    plt.scatter(X_axis + 0.2, rho_fitted_df.loc[target_choice_sets, 1], marker='^', color="green", facecolors='none', label=r"Fitted value ($x = 1$)", s=50)
    
    # plt.scatter(X_axis - 0.2, rho_df.loc[target_choice_sets, 0], marker='o', color="red", facecolors='none', label=r"True value ($x = 0$)", s=50)
    # plt.scatter(X_axis + 0.2, rho_df.loc[target_choice_sets, 1], marker='^', color="green", facecolors='none', label=r"True value ($x = 1$)", s=50)
    
    eb_0 = plt.errorbar(X_axis - 0.30, rum_bdd_mid_array_0, yerr=rum_error_array_0, capsize=5, fmt='o', markersize=0, ecolor='red', markeredgecolor="black", color='w', label=r"RUM bound ($x = 0$)")
    eb_1 = plt.errorbar(X_axis + 0.30, rum_bdd_mid_array_1, yerr=rum_error_array_1, capsize=5, fmt='o', markersize=0, ecolor='green', markeredgecolor="black", color='w', label=r"RUM bound ($x = 1$)")
    eb_naive = plt.errorbar(X_axis, naive_bdd_mid_array, yerr=naive_error_array, capsize=5, fmt='o', markersize=0, ecolor='blue', markeredgecolor="black", color='w', label="Naive bound", alpha=0.6)
    eb_mon =  plt.errorbar(X_axis - 0.15, mon_bdd_mid_array_0, yerr=mon_error_array_0, capsize=5, fmt='o', markersize=0, ecolor='purple', markeredgecolor="black", color='w', label=r"Monotonicity bound ($x = 0$)")
    eb_mon1 =  plt.errorbar(X_axis + 0.15, mon_bdd_mid_array_1, yerr=mon_error_array_1, capsize=5, fmt='o', markersize=0, ecolor='orange', markeredgecolor="black", color='w', label=r"Monotonicty bound ($x = 1$)")
    
    eb_0[-1][0].set_linestyle("solid")
    eb_mon[-1][0].set_linestyle("dashdot")
    eb_mon1[-1][0].set_linestyle("dotted")
    eb_1[-1][0].set_linestyle("dashed")
    eb_naive[-1][0].set_linestyle("dotted")

    plt.xticks(X_axis, ["{" + ", ".join(map(str, D)) + "}" for D in target_choice_sets])
    plt.xlabel("Choice set", fontsize=15)
    plt.ylabel("Choice Probability", fontsize=15)
    plt.legend(fontsize=12)
    plt.savefig("outputs/id_comparison.png", dpi=1000)
    # plt.show()