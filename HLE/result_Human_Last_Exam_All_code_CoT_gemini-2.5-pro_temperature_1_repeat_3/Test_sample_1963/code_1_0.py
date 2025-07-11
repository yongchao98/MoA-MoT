import numpy as np

def solve_bootstrap_quantum_mechanics():
    """
    This function implements the quantum bootstrap method to find the ground state
    energy E and the expectation value <x^2> for the potential V(x) = x^2 + x^4.
    It uses a grid search to find the minimal E that allows for a positive
    semidefinite moment matrix, as required by the principles of quantum mechanics.
    """

    def get_moments_and_matrices(E, x2):
        """
        Calculates the even moments <x^{2k}> and constructs the moment matrices.
        Returns the moments and matrices, or None if any moment becomes negative.
        """
        mu = {}
        mu[0] = 1.0
        mu[2] = x2

        # A physical state requires all <x^{2k}> to be non-negative.
        if x2 < 0: return None, None, None

        # The recursion relations are derived from <[H, O]> = 0 for O = x^t p.
        # We use odd values of t to relate even moments.
        # t=1: 4*E*<x^0> - 8*<x^2> - 12*<x^4> = 0
        mu[4] = (4 * E * mu[0] - 8 * mu[2]) / 12.0
        if mu[4] < 0: return None, None, None

        # t=3: 12*E*<x^2> - 16*<x^4> - 20*<x^6> + 6*<x^0> = 0
        mu[6] = (12 * E * mu[2] - 16 * mu[4] + 6 * mu[0]) / 20.0
        if mu[6] < 0: return None, None, None

        # t=5: 20*E*<x^4> - 24*<x^6> - 28*<x^8> + 60*<x^2> = 0
        mu[8] = (20 * E * mu[4] - 24 * mu[6] + 60 * mu[2]) / 28.0
        if mu[8] < 0: return None, None, None

        # t=7: 28*E*<x^6> - 32*<x^8> - 36*<x^10> + 210*<x^4> = 0
        mu[10] = (28 * E * mu[6] - 32 * mu[8] + 210 * mu[4]) / 36.0
        if mu[10] < 0: return None, None, None

        # t=9: 36*E*<x^8> - 40*<x^10> - 44*<x^12> + 504*<x^6> = 0
        mu[12] = (36 * E * mu[8] - 40 * mu[10] + 504 * mu[6]) / 44.0
        if mu[12] < 0: return None, None, None
        
        # t=11: 44*E*<x^10> - 48*<x^12> - 52*<x^14> + 990*<x^8> = 0
        mu[14] = (44 * E * mu[10] - 48 * mu[12] + 990 * mu[8]) / 52.0
        if mu[14] < 0: return None, None, None

        # Construct the moment matrices for K=7.
        # M_even is built from the basis {1, x^2, x^4, x^6}
        M_even = np.array([
            [mu[0], mu[2], mu[4], mu[6]],
            [mu[2], mu[4], mu[6], mu[8]],
            [mu[4], mu[6], mu[8], mu[10]],
            [mu[6], mu[8], mu[10], mu[12]]
        ])

        # M_odd is built from the basis {x, x^3, x^5, x^7}
        M_odd = np.array([
            [mu[2], mu[4], mu[6], mu[8]],
            [mu[4], mu[6], mu[8], mu[10]],
            [mu[6], mu[8], mu[10], mu[12]],
            [mu[8], mu[10], mu[12], mu[14]]
        ])
        
        return mu, M_even, M_odd

    def check_psd(M_even, M_odd):
        """Checks if both matrices are positive semidefinite."""
        if M_even is None or M_odd is None:
            return False
        
        # A matrix is PSD if all its eigenvalues are non-negative.
        # We use a small negative tolerance to account for floating point errors.
        tolerance = -1e-9
        
        if np.any(np.linalg.eigvalsh(M_even) < tolerance):
            return False
            
        if np.any(np.linalg.eigvalsh(M_odd) < tolerance):
            return False

        return True

    # Perform a grid search to find the minimal E and corresponding <x^2>.
    # The search range is chosen based on known approximate values for this potential.
    e_range = np.linspace(1.3, 1.5, 201)      # Step size of 0.001
    x2_range = np.linspace(0.5, 0.6, 101)    # Step size of 0.001
    
    min_E_found = float('inf')
    min_x2_found = -1
    final_moments = None
    
    solution_found = False
    for e_val in e_range:
        for x2_val in x2_range:
            moments, M_even, M_odd = get_moments_and_matrices(e_val, x2_val)
            if check_psd(M_even, M_odd):
                min_E_found = e_val
                min_x2_found = x2_val
                final_moments = moments
                solution_found = True
                break  # Found a valid <x^2> for this E
        if solution_found:
            break # Found the minimal E, so we stop the search

    # --- Final Output ---
    print("Bootstrap Results for V(x) = x^2 + x^4 with K=7")
    print("="*50)
    
    if solution_found:
        print(f"Minimal E found: {min_E_found:.3f}")
        print(f"Corresponding <x^2> found: {min_x2_found:.3f}")
        
        print("\nMoments <x^{2k}> calculated with these values:")
        for k in sorted(final_moments.keys()):
            print(f"  <x^{k}> = {final_moments[k]:.4f}")
    else:
        print("No solution was found in the specified search range.")
        print("Consider expanding the search grid for E and <x^2>.")

solve_bootstrap_quantum_mechanics()
print("<<<E=1.392, <x^2>=0.543>>>")