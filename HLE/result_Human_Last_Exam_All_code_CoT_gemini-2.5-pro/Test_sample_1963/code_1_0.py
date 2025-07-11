import numpy as np

def solve_bootstrap():
    """
    Implements the quantum bootstrap method to find the ground state energy E
    and the expectation value <x^2> for the potential V(x) = x^2 + x^4.
    """

    def calculate_moments(E, x2, K):
        """
        Calculates even moments <x^{2n}> up to 2K using the recursion relation.
        Returns a dictionary of moments or None if a physical constraint is violated.
        """
        moments = {}
        # Due to V(x) being even, all odd moments <x^{2n+1}> are 0.
        # We start with the initial conditions.
        moments[0] = 1.0
        moments[2] = x2

        # A key physical constraint is that <x^{2n}> must be non-negative.
        # The first recursion gives <x^4> = (E - 2*<x^2>)/3.
        # If <x^4> is negative, the state is unphysical. This prunes the search space.
        if (E - 2 * x2) < 0:
            return None
        moments[4] = (E - 2 * x2) / 3.0

        # Now, we use the general recursion relation derived from the bootstrap equations:
        # (4t+8)<x^{t+3}> = 4tE<x^{t-1}> + t(t-1)(t-2)<x^{t-3}> - (4t+4)<x^{t+1}>
        # To get relations between even moments, we must use odd values for t.
        # We need moments up to <x^{2K}> = <x^{14}>.
        
        # Loop for t = 3, 5, 7, 9, 11
        for i in range(1, K - 1):
            t = float(2 * i + 1)
            
            new_moment_idx = int(t + 3) # e.g., t=3 -> idx=6
            
            den = (4 * t + 8)
            
            term1 = 4 * t * E * moments[int(t - 1)]
            term2 = t * (t - 1) * (t - 2) * moments[int(t - 3)]
            term3 = - (4 * t + 4) * moments[int(t + 1)]
            
            moments[new_moment_idx] = (term1 + term2 + term3) / den
            
            # Check for positivity of moments during calculation
            if moments[new_moment_idx] < 0:
                return None

        return moments

    def check_positive_semidefinite(E, x2, K=7):
        """
        Calculates moments and checks if the corresponding matrices are positive semidefinite.
        """
        moments = calculate_moments(E, x2, K)
        if moments is None:
            return False

        # For K=7, the basis for operators splits into even and odd polynomials.
        # Even basis: {1, x^2, x^4, x^6}, giving a 4x4 matrix M_e.
        # Odd basis: {x, x^3, x^5, x^7}, giving a 4x4 matrix M_o.

        try:
            # M_e[i,j] = <x^{2i+2j}> for i,j in {0,1,2,3}
            M_e = np.array([
                [moments[0], moments[2], moments[4], moments[6]],
                [moments[2], moments[4], moments[6], moments[8]],
                [moments[4], moments[6], moments[8], moments[10]],
                [moments[6], moments[8], moments[10], moments[12]]
            ])
            
            # M_o[i,j] = <x^{2i+2j+2}> for i,j in {0,1,2,3}
            M_o = np.array([
                [moments[2], moments[4], moments[6], moments[8]],
                [moments[4], moments[6], moments[8], moments[10]],
                [moments[6], moments[8], moments[10], moments[12]],
                [moments[8], moments[10], moments[12], moments[14]]
            ])
        except KeyError:
            # A required moment was not calculated.
            return False

        # A matrix is positive semidefinite if all its eigenvalues are non-negative.
        # We use a small tolerance for floating point inaccuracies.
        epsilon = -1e-9
        if np.any(np.linalg.eigvalsh(M_e) < epsilon):
            return False
        if np.any(np.linalg.eigvalsh(M_o) < epsilon):
            return False
            
        return True

    # Perform a grid search to find the minimal allowed (E, <x^2>).
    # We search a reasonable parameter space with a fine step.
    min_E_found = None
    min_x2_found = None

    E_range = np.arange(1.3, 1.5, 0.001)
    x2_range = np.arange(0.3, 0.5, 0.001)

    for E in E_range:
        for x2 in x2_range:
            if check_positive_semidefinite(E, x2):
                min_E_found = E
                min_x2_found = x2
                # Since we are iterating from low to high, the first valid point
                # found is the minimum on our grid.
                break
        if min_E_found is not None:
            break

    if min_E_found is not None and min_x2_found is not None:
        print("Found minimal values based on the bootstrap constraints:")
        print(f"Minimal E = {min_E_found:.3f}")
        print(f"Minimal <x^2> = {min_x2_found:.3f}")
    else:
        print("No solution found in the specified range.")
        print("Please consider expanding the search ranges for E and <x^2>.")

solve_bootstrap()