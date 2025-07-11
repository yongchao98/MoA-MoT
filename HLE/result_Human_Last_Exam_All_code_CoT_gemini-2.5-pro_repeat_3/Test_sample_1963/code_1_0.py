import numpy as np

def find_ground_state_bootstrap():
    """
    Performs a grid search to find the minimal E and corresponding <x^2>
    for the potential V(x) = x^2 + x^4 using the quantum bootstrap method with K=7.
    """

    def check_psd(E, X2):
        """
        Calculates moments and checks if the Hankel matrices are positive semidefinite.
        Returns True if the point (E, X2) is allowed, False otherwise.
        """
        # The moments dict will store X_n = <x^n>
        X = {}
        X[0] = 1.0
        X[2] = X2

        # Check for non-physical input
        if X2 < 0:
            return False

        # Recursion relations derived from the bootstrap equations for V(x) = x^2 + x^4
        # t=1 -> X4 = (E - 2*X2) / 3
        X[4] = (E - 2 * X[2]) / 3.0
        
        # t=3 -> X6 = (12*E*X2 - 16*X4 + 6) / 20
        X[6] = (12 * E * X[2] - 16 * X[4] + 6) / 20.0
        
        # t=5 -> X8 = (5*E*X4 - 6*X6 + 15*X2) / 7
        X[8] = (5 * E * X[4] - 6 * X[6] + 15 * X[2]) / 7.0
        
        # t=7 -> X10 = (14*E*X6 - 16*X8 + 105*X4) / 18
        X[10] = (14 * E * X[6] - 16 * X[8] + 105 * X[4]) / 18.0
        
        # t=9 -> X12 = (9*E*X8 - 10*X10 + 126*X6) / 11
        X[12] = (9 * E * X[8] - 10 * X[10] + 126 * X[6]) / 11.0
        
        # t=11 -> X14 = (22*E*X10 - 24*X12 + 495*X8) / 26
        X[14] = (22 * E * X[10] - 24 * X[12] + 495 * X[8]) / 26.0

        # All even moments must be non-negative
        if any(val < -1e-9 for val in X.values()):
             return False

        # Construct the two Hankel matrices for K=7
        # M_even corresponds to the basis {1, x^2, x^4, x^6}
        M_even = np.array([
            [X[0], X[2], X[4], X[6]],
            [X[2], X[4], X[6], X[8]],
            [X[4], X[6], X[8], X[10]],
            [X[6], X[8], X[10], X[12]]
        ])

        # M_odd corresponds to the basis {x, x^3, x^5, x^7}
        M_odd = np.array([
            [X[2], X[4], X[6], X[8]],
            [X[4], X[6], X[8], X[10]],
            [X[6], X[8], X[10], X[12]],
            [X[8], X[10], X[12], X[14]]
        ])

        try:
            # Check for positive semidefiniteness by ensuring all eigenvalues are non-negative.
            # A small tolerance is used for floating-point inaccuracies.
            tolerance = -1e-9
            eig_even = np.linalg.eigvalsh(M_even)
            if np.any(eig_even < tolerance):
                return False
            
            eig_odd = np.linalg.eigvalsh(M_odd)
            if np.any(eig_odd < tolerance):
                return False

        except np.linalg.LinAlgError:
            return False
            
        return True

    # Define the search grid. Ranges are chosen around known approximate values.
    # The grid resolution is chosen to provide at least 3-digit precision.
    e_range = np.linspace(1.38, 1.41, 151)
    x2_range = np.linspace(0.37, 0.39, 101)
    
    allowed_points = []
    
    print("Searching for allowed (E, <x^2>) pairs...")
    for E_test in e_range:
        for x2_test in x2_range:
            if check_psd(E_test, x2_test):
                allowed_points.append((E_test, x2_test))

    if not allowed_points:
        print("No allowed points were found in the specified search grid.")
        print("Please consider expanding the search ranges for E and <x^2>.")
        return

    # Find the point with the minimum energy from the list of allowed points
    min_E, best_x2 = min(allowed_points, key=lambda p: p[0])

    print("\n--- Bootstrap Calculation Results ---")
    print(f"The minimal ground state energy E found is: {min_E:.3f}")
    print(f"The corresponding expectation value <x^2> is: {best_x2:.3f}")
    
    # The final answer format as requested.
    # Using f-string for formatted output within the tag.
    final_answer = f"E={min_E:.3f}, <x^2>={best_x2:.3f}"
    print(f"\n<<<{final_answer}>>>")


if __name__ == '__main__':
    find_ground_state_bootstrap()
