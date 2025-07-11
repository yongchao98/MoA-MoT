import numpy as np

def solve_bootstrap_qm():
    """
    Solves for the ground state energy E and <x^2> of the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """

    def calculate_moments(E, m2):
        """
        Calculates even moments <x^{2k}> from k=0 to 7 using the recursion relation.
        
        Args:
            E (float): The test energy value.
            m2 (float): The test value for <x^2>.

        Returns:
            numpy.ndarray or None: An array of moments m[i] = <x^i>, or None if invalid.
        """
        # We need even moments m_0, m_2, ..., m_14 for K=7
        # The highest index needed is for M_o, which is <x^{14}>
        m = np.zeros(16)
        
        # Physical constraint: <x^2> must be non-negative
        if m2 < 0:
            return None

        # Initial conditions
        m[0] = 1.0  # <x^0> = 1
        m[2] = m2   # <x^2>

        # The recursion relation for the potential V(x) = x^2 + x^4 is:
        # (8k+12)m_{2k+4} = (8k+4)E*m_{2k} - (8k+8)m_{2k+2} + (2k+1)(2k)(2k-1)m_{2k-2}
        
        # Calculate m[4] (for k=0)
        # 12*m[4] = 4*E*m[0] - 8*m[2]
        m[4] = (4 * E * m[0] - 8 * m[2]) / 12.0
        
        # Calculate higher moments using a loop for k = 1, 2, ..., 6
        for k in range(1, 7):
            idx_2k_minus_2 = 2 * k - 2
            idx_2k = 2 * k
            idx_2k_plus_2 = 2 * k + 2
            idx_2k_plus_4 = 2 * k + 4
            
            term_E = (8 * k + 4) * E * m[idx_2k]
            term_m_plus_2 = (8 * k + 8) * m[idx_2k_plus_2]
            term_m_minus_2 = (2 * k + 1) * (2 * k) * (2 * k - 1) * m[idx_2k_minus_2]
            
            numerator = term_E - term_m_plus_2 + term_m_minus_2
            denominator = 8 * k + 12
            
            m[idx_2k_plus_4] = numerator / denominator
            
        return m

    def check_psd(moments):
        """
        Checks if the moment matrices M_e and M_o are positive semidefinite.

        Args:
            moments (numpy.ndarray): Array of moments.

        Returns:
            bool: True if both matrices are positive semidefinite, False otherwise.
        """
        if moments is None:
            return False
        
        m = moments
        
        # Construct the even-power matrix M_e (4x4 for K=7)
        # M_e[i,j] = <x^{2(i+j)}> for i,j in {0,1,2,3}
        M_e = np.array([
            [m[0], m[2], m[4], m[6]],
            [m[2], m[4], m[6], m[8]],
            [m[4], m[6], m[8], m[10]],
            [m[6], m[8], m[10], m[12]]
        ])
        
        # Construct the odd-power matrix M_o (4x4 for K=7)
        # M_o[i,j] = <x^{2(i+j+1)}> for i,j in {0,1,2,3}
        M_o = np.array([
            [m[2], m[4], m[6], m[8]],
            [m[4], m[6], m[8], m[10]],
            [m[6], m[8], m[10], m[12]],
            [m[8], m[10], m[12], m[14]]
        ])
        
        # A symmetric matrix is positive semidefinite if all its eigenvalues are non-negative.
        # We use a small tolerance to account for floating point inaccuracies.
        tolerance = -1e-9
        
        try:
            eig_e = np.linalg.eigvalsh(M_e)
            eig_o = np.linalg.eigvalsh(M_o)
        except np.linalg.LinAlgError:
            # Matrix may not be valid (e.g., contains NaNs or Infs)
            return False
            
        return np.all(eig_e >= tolerance) and np.all(eig_o >= tolerance)

    # --- Main Search Logic ---
    min_E_found = float('inf')
    best_m2_for_min_E = float('inf')

    # Define a reasonable search grid based on physical estimates (e.g., perturbation theory)
    E_range = np.arange(1.3, 1.5, 0.001)
    m2_range = np.arange(0.3, 0.5, 0.001)

    # Iterate through E values from low to high
    for E_test in E_range:
        found_valid_m2 = False
        
        # For a given E, find the minimum valid m2
        for m2_test in m2_range:
            moments = calculate_moments(E_test, m2_test)
            if check_psd(moments):
                # We found a valid (E, m2) pair.
                # Since we are iterating E from low to high, this E_test is the minimal energy.
                # The first valid m2_test will be the minimal one for this E.
                min_E_found = E_test
                best_m2_for_min_E = m2_test
                found_valid_m2 = True
                break # Found the minimal m2 for this E, move to printing results
        
        if found_valid_m2:
            break # Found the minimal E, exit the main search loop

    # --- Output the final results ---
    print("Found minimal values based on the bootstrap method with K=7.")
    print(f"Minimal E = {min_E_found:.3f}")
    print(f"Minimal <x^2> = {best_m2_for_min_E:.3f}")

# Execute the function
solve_bootstrap_qm()