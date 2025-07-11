import math

def nCr_q(n, r, q):
    """Calculates the q-binomial coefficient (Gaussian coefficient)."""
    if r < 0 or r > n:
        return 0
    num = 1
    den = 1
    for i in range(r):
        num *= (q**(n - i) - 1)
        den *= (q**(i + 1) - 1)
    return num // den

def GL_size(n, q):
    """Calculates the order of the general linear group GL_n(q)."""
    if n == 0:
        return 1
    size = 1
    for i in range(n):
        size *= (q**n - q**i)
    return size

def calculate_proportion():
    """
    Calculates the proportion of irreducible (3, 2)-stingray duos in GL_5(4) x GL_5(4).
    """
    q = 4
    d = 5
    e1 = 3
    e2 = 2
    
    dim_hom_U2_U1 = e2 * e1 # dim Hom(U2, U1)
    dim_hom_U1_U2 = e1 * e2 # dim Hom(U1, U2)

    # Number of linear maps alpha: U2 -> U1
    num_alpha = q**dim_hom_U2_U1
    # Number of linear maps beta: U1 -> U2
    num_beta = q**dim_hom_U1_U2

    # --- Count maps alpha by rank ---
    # Number of maps alpha of rank 1
    # N(k,m,r) = [m,r]_q * Product_{i=0..r-1} (q^k - q^i)
    # k=e2=2, m=e1=3, r=1
    num_alpha_r1 = nCr_q(e1, 1, q) * (q**e2 - 1)

    # Number of maps alpha of rank 2
    # k=e2=2, m=e1=3, r=2
    num_alpha_r2 = nCr_q(e1, 2, q) * (q**e2 - 1) * (q**e2 - q)
    
    # Check
    assert num_alpha_r1 + num_alpha_r2 == num_alpha - 1

    # --- Count "bad" pairs (alpha, beta) for reducibility condition (3) ---
    # This is the set Z where det(I - beta*alpha) = 0
    
    # Case alpha has rank 1:
    # beta*alpha has rank at most 1. det(I - M) = 0 iff tr(M)=1 for rank 1 M.
    # The map beta -> tr(beta*alpha) is a non-zero linear functional.
    # Its kernel has size q^(dim-1). The number of solutions to tr=1 is also q^(dim-1).
    num_beta_bad_for_alpha_r1 = q**(dim_hom_U1_U2 - 1)
    
    # Case alpha has rank 2:
    # beta*alpha can be any 2x2 matrix with appropriate choice of beta.
    # The number of beta for each such 2x2 matrix M is q^(e1*e2 - e2*e2).
    # Number of singular 2x2 matrices is q^(2*2) - |GL(2,q)|.
    # The number of 2x2 matrices M with det(I-M)=0 is same as singular matrices.
    num_singular_M2 = q**(e2*e2) - GL_size(e2, q)
    num_beta_factor = q**(e1*e2 - e2*e2)
    num_beta_bad_for_alpha_r2 = num_singular_M2 * num_beta_factor

    # Total size of Z (pairs (alpha!=0, beta) with det=0)
    # The condition det(I - beta*alpha)=0 implies beta!=0
    num_Z = (num_alpha_r1 * num_beta_bad_for_alpha_r1) + \
              (num_alpha_r2 * num_beta_bad_for_alpha_r2)
              
    # --- Total number of reducible configurations ---
    # We established the three conditions are mutually exclusive if we formulate them carefully.
    # 1. alpha = 0: 1 * num_beta = q^6 pairs
    # 2. beta = 0 (and alpha != 0): (num_alpha - 1) * 1 = q^6 - 1 pairs
    # 3. det(I - beta*alpha) = 0 (and alpha != 0, beta != 0): num_Z pairs
    # The logic used for num_Z already assumed alpha!=0 and showed beta!=0.
    
    num_bad_configs_alpha_eq_0 = num_beta
    num_bad_configs_beta_eq_0 = num_alpha - 1
    num_bad_configs_det_eq_0 = num_Z

    num_bad_configs = num_bad_configs_alpha_eq_0 + num_bad_configs_beta_eq_0 + num_bad_configs_det_eq_0
    
    total_configs = num_alpha * num_beta
    num_good_configs = total_configs - num_bad_configs

    proportion = num_good_configs / total_configs
    
    print("(a) No")
    print("(b) {(1), (2), (3)}")
    print("(c) The proportion is calculated as (Total Configurations - Reducible Configurations) / Total Configurations.")
    print(f"Total Configurations = {total_configs}")
    print(f"Reducible due to U2 = F1 (alpha=0): {num_bad_configs_alpha_eq_0}")
    print(f"Reducible due to U1 = F2 (beta=0, alpha!=0): {num_bad_configs_beta_eq_0}")
    print(f"Reducible due to F1_intersect_F2!=0 (alpha,beta!=0): {num_bad_configs_det_eq_0}")
    print(f"Total Reducible Configurations = {num_bad_configs}")
    print(f"Irreducible Configurations = {num_good_configs}")
    print(f"Proportion = {num_good_configs} / {total_configs}")
    print(f"Final Proportion Value = {proportion}")

calculate_proportion()
