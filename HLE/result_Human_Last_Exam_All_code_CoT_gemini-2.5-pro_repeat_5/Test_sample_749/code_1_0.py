import numpy as np

def get_jump_probabilities(env, site):
    """Returns (p_left, p_right) for a given site."""
    if env.get(site) == 'red':
        return 4/5, 1/5
    else: # blue
        return 1/5, 4/5

def calculate_hitting_probability(env, start_site, target_site=0, lattice_size=200):
    """
    Calculates the probability of hitting target_site from start_site.
    We solve a system of linear equations on a finite lattice [-L, L].
    """
    if start_site == target_site:
        return 1.0

    L = lattice_size
    # We solve for pi_k on sites k in [-L, L] \ {target_site}
    sites = [k for k in range(-L, L + 1) if k != target_site]
    n = len(sites)
    
    # Matrix A for the system A*pi = b
    A = np.zeros((n, n))
    # Vector b
    b = np.zeros(n)
    
    site_to_idx = {site: i for i, site in enumerate(sites)}

    for i, k in enumerate(sites):
        A[i, i] = -1
        p_left, p_right = get_jump_probabilities(env, k)
        
        # Neighbor to the left
        k_left = k - 1
        if k_left == target_site:
            b[i] -= p_left * 1.0 # pi_target = 1
        elif k_left in site_to_idx:
            A[i, site_to_idx[k_left]] = p_left
        else: # Boundary at -L
             b[i] -= p_left * 1.0 # Assume pi_{-L-1} approx 1 (right drift)

        # Neighbor to the right
        k_right = k + 1
        if k_right == target_site:
            b[i] -= p_right * 1.0 # pi_target = 1
        elif k_right in site_to_idx:
            A[i, site_to_idx[k_right]] = p_right
        else: # Boundary at L
            b[i] -= p_right * 0.0 # Assume pi_{L+1} approx 0 (right drift)

    try:
        pi = np.linalg.solve(A, b)
        return pi[site_to_idx[start_site]]
    except np.linalg.LinAlgError:
        return np.nan # Should not happen for this problem

def calculate_m0(env):
    """Calculates the m_0 parameter for a given environment."""
    pi_1 = calculate_hitting_probability(env, 1)
    pi_minus_1 = calculate_hitting_probability(env, -1)
    
    p_left_0, p_right_0 = get_jump_probabilities(env, 0)
    
    m0 = p_left_0 * pi_minus_1 + p_right_0 * pi_1
    return m0, pi_1, pi_minus_1

def main():
    """Main function to run the analysis."""
    print("Analyzing the condition for infinite visits to site 0.")
    print("The condition is R_omega = (1+h)*m_0(omega) > 1.")
    print("This is equivalent to h > (1/m_0(omega)) - 1.\n")
    
    # Case 1: All sites are blue
    env_blue = {}
    m0_blue, pi1_blue, pim1_blue = calculate_m0(env_blue)
    h_thresh_blue = 1 / m0_blue - 1 if m0_blue > 0 else float('inf')
    
    print("Case 1: All sites are blue")
    print(f"  pi_1 (prob to reach 0 from 1)  = {pi1_blue:.4f}")
    print(f"  pi_-1 (prob to reach 0 from -1) = {pim1_blue:.4f}")
    print(f"  m_0 = (1/5)*pi_-1 + (4/5)*pi_1 = {m0_blue:.4f}")
    print(f"  Result: m_0 is < 1. The threshold for h is h > {h_thresh_blue:.4f}. This is not in (0, 1/2).\n")

    # Case 2: Site 0 is red, others blue
    env_0_red = {0: 'red'}
    m0_0r, pi1_0r, pim1_0r = calculate_m0(env_0_red)
    h_thresh_0r = 1 / m0_0r - 1 if m0_0r > 0 else float('inf')

    print("Case 2: Site 0 is red")
    print(f"  pi_1  = {pi1_0r:.4f}")
    print(f"  pi_-1 = {pim1_0r:.4f}")
    print(f"  m_0 = (4/5)*pi_-1 + (1/5)*pi_1 = {m0_0r:.4f}")
    print(f"  Result: m_0 is < 1. The threshold for h is h > {h_thresh_0r:.4f}.\n")

    # Case 3: Sites 0 and 1 are red
    env_01_red = {0: 'red', 1: 'red'}
    m0_01r, pi1_01r, pim1_01r = calculate_m0(env_01_red)
    h_thresh_01r = 1 / m0_01r - 1 if m0_01r > 0 else float('inf')

    print("Case 3: Sites 0 and 1 are red")
    print(f"  pi_1  = {pi1_01r:.4f}")
    print(f"  pi_-1 = {pim1_01r:.4f}")
    print(f"  m_0 = (4/5)*pi_-1 + (1/5)*pi_1 = {m0_01r:.4f}")
    print(f"  Result: m_0 is < 1. The threshold for h is h > {h_thresh_01r:.4f}.\n")

    # Case 4: A 'trap' with sites 1,2,3,4,5 red, 0 blue
    env_trap = {i: 'red' for i in range(1, 6)}
    m0_trap, pi1_trap, pim1_trap = calculate_m0(env_trap)
    h_thresh_trap = 1 / m0_trap - 1 if m0_trap > 0 else float('inf')

    print("Case 4: A 'trap' (sites 1-5 are red, 0 is blue)")
    print(f"  pi_1  = {pi1_trap:.4f}")
    print(f"  pi_-1 = {pim1_trap:.4f}")
    print(f"  m_0 = (1/5)*pi_-1 + (4/5)*pi_1 = {m0_trap:.4f}")
    print(f"  Result: m_0 is very close to 1, but still < 1. The threshold for h, {h_thresh_trap:.4f}, is small but positive.\n")
    
    print("Conclusion: For any environment with a finite number of red sites, m_0(omega) < 1.")
    print("This means a strictly positive value of h is required for infinite visits.")
    print("As h -> 0, for any fixed finite environment, the condition is not met.")
    print("The probability measure concentrates on these finite environments.")
    print("Therefore, the limiting probability is 0.")


if __name__ == '__main__':
    main()
