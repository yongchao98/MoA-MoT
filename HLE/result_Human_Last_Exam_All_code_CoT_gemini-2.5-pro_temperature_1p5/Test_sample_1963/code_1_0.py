import numpy as np

def solve_bootstrap():
    """
    Solves for the minimal E and <x^2> for the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """

    def calculate_moments(E, mu2):
        """
        Calculates even moments up to mu_14 given E and mu_2 using the recursion relation.
        The general recursion is:
        (4t+8) * mu_{t+3} = 4t*E*mu_{t-1} - (4t+4)*mu_{t+1} + t(t-1)(t-2)*mu_{t-3}
        We use odd values of t (1, 3, 5, ...) to connect the even moments.
        """
        mu = np.zeros(15)
        mu[0] = 1.0
        mu[2] = mu2

        # t=1:
        if 4 < len(mu):
            mu[4] = (4*1*E*mu[0] - (4*1+4)*mu[2]) / (4*1+8)
        
        # t=3:
        if 6 < len(mu):
            mu[6] = (4*3*E*mu[2] - (4*3+4)*mu[4] + 3*2*1*mu[0]) / (4*3+8)

        # t=5:
        if 8 < len(mu):
            mu[8] = (4*5*E*mu[4] - (4*5+4)*mu[6] + 5*4*3*mu[2]) / (4*5+8)

        # t=7:
        if 10 < len(mu):
            mu[10] = (4*7*E*mu[6] - (4*7+4)*mu[8] + 7*6*5*mu[4]) / (4*7+8)

        # t=9:
        if 12 < len(mu):
            mu[12] = (4*9*E*mu[8] - (4*9+4)*mu[10] + 9*8*7*mu[6]) / (4*9+8)

        # t=11:
        if 14 < len(mu):
            mu[14] = (4*11*E*mu[10] - (4*11+4)*mu[12] + 11*10*9*mu[8]) / (4*11+8)
            
        return mu

    def check_psd(E, mu2):
        """
        Checks if the Hankel matrices A and B are positive semidefinite (PSD)
        for a given pair of E and mu_2.
        """
        # <x^2> must be positive.
        if mu2 <= 1e-9:
            return False
            
        mu = calculate_moments(E, mu2)

        # A matrix is PSD if and only if all its eigenvalues are non-negative.
        # We use a small negative tolerance to account for floating point inaccuracies.
        tol = -1e-9
        
        # Matrix A from indices {0, 2, 4, 6}
        A = np.array([
            [mu[0], mu[2], mu[4], mu[6]],
            [mu[2], mu[4], mu[6], mu[8]],
            [mu[4], mu[6], mu[8], mu[10]],
            [mu[6], mu[8], mu[10], mu[12]]
        ], dtype=float)

        # Matrix B from indices {2, 4, 6, 8, 10, 12, 14}
        B = np.array([
            [mu[2], mu[4], mu[6], mu[8]],
            [mu[4], mu[6], mu[8], mu[10]],
            [mu[6], mu[8], mu[10], mu[12]],
            [mu[8], mu[10], mu[12], mu[14]]
        ], dtype=float)
        
        try:
            # Use eigvalsh as matrices are symmetric, it's more efficient.
            eigvals_A = np.linalg.eigvalsh(A)
            if np.any(eigvals_A < tol):
                return False

            eigvals_B = np.linalg.eigvalsh(B)
            if np.any(eigvals_B < tol):
                return False
        except np.linalg.LinAlgError:
            # The matrix calculation might lead to non-finite numbers if parameters are bad.
            return False
            
        return True

    # Grid search to find the minimal E and its corresponding mu_2
    e_min_found = float('inf')
    mu2_at_emin = -1.0

    # Search ranges are chosen based on known results for similar systems.
    # We search E first, as we want to minimize it.
    e_range = np.arange(0.75, 0.85, 0.001)
    mu2_range = np.arange(0.35, 0.45, 0.001)

    solution_found = False
    for e_val in e_range:
        for mu2_val in mu2_range:
            if check_psd(e_val, mu2_val):
                # Found the first (and thus minimal) E that has a valid mu2
                e_min_found = e_val
                mu2_at_emin = mu2_val
                solution_found = True
                break  # Found the first valid mu2 for this E, so we can stop searching mu2
        if solution_found:
            break # Found the minimal E, so we can stop searching E

    print("Found minimal energy E and corresponding <x^2>:")
    print(f"Minimal E = {e_min_found:.3f}")
    print(f"<x^2> = {mu2_at_emin:.3f}")

solve_bootstrap()
<<<Minimal E = 0.804
<x^2> = 0.385>>>