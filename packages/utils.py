from itertools import chain, combinations


def get_power_sets(universal):
    """
    Parameters
    ----------
    universal : iterable
        universal set

    Returns
    -------
    list_power_sets : list
        list of all the power sets

    Examples
    --------
    >>> s = [1, 2, 3]
    >>> get_power_sets(s)
    [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    """
    s = list(universal)
    list_power_sets = list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
    return list_power_sets


def from_df_to_func(rho_df):
    return lambda D, x: rho_df.loc[[D], x][0]

