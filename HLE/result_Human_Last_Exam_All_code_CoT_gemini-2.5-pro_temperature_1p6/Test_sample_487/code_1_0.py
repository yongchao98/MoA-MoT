def calculate_h11_max():
    """
    This function calculates the maximal value of the Hodge number h^{1,1}(M)
    based on the analysis described above. The problem is a constrained integer
    optimization problem.

    Let n_c be the number of curve components of the fixed locus S^rho.
    Let g_tot be the total genus of these curves.
    Let m be the number of isolated fixed points.
    We found that m must be even, so let m = 2k.
    k_psi is the number of fixed points of the involution on C. We choose k_psi=6.

    We want to maximize h = 11 + 7*n_c - g_tot + 13*k
    Subject to constraints:
    1) n_c, g_tot, k >= 0 (integers)
    2) n_c + k <= g_tot  (from r <= 10 for disconnected locus)
    3) n_c + k + g_tot <= 12 (from Betti number bound)
    """

    max_h = 0
    best_config = {}

    # We are maximizing a function linear in n_c, k and -g_tot.
    # To maximize it, we should maximize k and n_c and minimize g_tot.
    # From constraint (2), g_tot >= n_c + k. Let's start by trying to find the
    # maximum possible k, which has the largest coefficient.

    # From (2) and (3):
    # (n_c + k) + g_tot <= 12
    # Since g_tot >= n_c + k, we have
    # (n_c + k) + (n_c + k) <= 12  =>  2*n_c + 2*k <= 12  => n_c + k <= 6
    # To maximize k, we should minimize n_c.
    # A fixed locus with curves requires n_c >= 1.
    
    n_c = 1
    # Now, with n_c=1, the constraint becomes 1 + k <= 6, so k <= 5.
    # Let's take the maximum possible value for k.
    k = 5
    
    # Now find g_tot. From constraint (2), at its boundary:
    # g_tot = n_c + k
    g_tot = n_c + k
    
    # Check if this satisfies constraint (3):
    # n_c + k + g_tot <= 12
    # 1 + 5 + 6 = 12, which satisfies the constraint.
    
    # Now we have found the optimal configuration.
    # n_c = 1
    # k = 5, which means m = 2*k = 10 isolated points
    # g_tot = 6, for the single curve component
    
    config = {'n_c': n_c, 'g_tot': g_tot, 'k_m': k, 'm': 2*k}

    # Now we calculate the Hodge number h^{1,1}(M) for this configuration.
    h_1_1_M = 11 + 7 * config['n_c'] - config['g_tot'] + 13 * config['k_m']
    
    # For clarity, let's also calculate it with the original formula.
    # r = 10 + n_c - g_tot + k
    r = 10 + config['n_c'] - config['g_tot'] + config['k_m']
    # c_rho = n_c + m
    c_rho = config['n_c'] + config['m']
    k_psi = 6
    h_1_1_M_alt = r + 1 + k_psi * c_rho
    
    print("Optimal configuration for the fixed locus S^rho:")
    print(f"  Number of curve components (n_c): {config['n_c']}")
    print(f"  Genus of the curve (g_tot): {config['g_tot']}")
    print(f"  Number of isolated points (m): {config['m']}")
    
    print("\nDerived parameters:")
    print(f"  Rank of invariant lattice on H^1,1 (r): {r}")
    print(f"  Number of components of S^rho (c_rho): {c_rho}")
    print(f"  Number of fixed points on C (k_psi): {k_psi}")

    print("\nCalculation of the maximal Hodge number h^{1,1}(M):")
    print(f"h^{{1,1}}(M) = r + 1 + k_psi * c_rho")
    print(f"           = {r} + 1 + {k_psi} * {c_rho}")
    print(f"           = {r+1} + {k_psi*c_rho}")
    print(f"           = {h_1_1_M_alt}")
    
calculate_h11_max()
