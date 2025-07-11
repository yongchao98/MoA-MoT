import numpy as np

def solve_bootstrap_qm():
    """
    Finds the minimal E and <x^2> for V(x) = x^2 + x^4 using the bootstrap method.
    """

    def calculate_moments(E, m1):
        """
        Calculates even moments m_k = <x^{2k}> for k=0 to 7 based on E and m1.
        
        The recursion relation is derived from the bootstrap equations for V(x) = x^2 + x^4:
        4(2k+3) * m_{k+2} = (2k+1)(2k)(2k-1) * m_{k-1} + 4(2k+1)E * m_k - 4(2k+2) * m_{k+1}
        """
        if m1 < 0:
            return None
            
        m = np.zeros(8)
        m[0] = 1.0  # <x^0> = 1
        m[1] = m1    # <x^2>
        
        try:
            # For t=1 (k=0 in our formula notation): E - 2*m1 - 3*m2 = 0
            m[2] = (E - 2 * m[1]) / 3.0

            # k=1 for m_3 = <x^6>
            m[3] = (3 * E * m[1] - 4 * m[2] + 1.5 * m[0]) / 5.0
            
            # k=2 for m_4 = <x^8>
            m[4] = (5 * E * m[2] - 6 * m[3] + 15 * m[1]) / 7.0
            
            # k=3 for m_5 = <x^10>
            m[5] = (7 * E * m[3] - 8 * m[4] + 52.5 * m[2]) / 9.0

            # k=4 for m_6 = <x^12>
            m[6] = (9 * E * m[4] - 10 * m[5] + 126 * m[3]) / 11.0

            # k=5 for m_7 = <x^14>
            m[7] = (11 * E * m[5] - 12 * m[6] + 247.5 * m[4]) / 13.0
        except (OverflowError, ValueError):
            return None # Invalid numbers encountered
            
        return m

    def check_psd(E, m1):
        """
        Checks if the Hankel matrices M_even and M_odd are positive semidefinite.
        """
        moments = calculate_moments(E, m1)
        if moments is None or np.any(np.isnan(moments)) or np.any(np.isinf(moments)):
            return False

        m = moments
        
        # M_ij = <x^{i+j}>, i,j from {0,2,4,6}
        M_even = np.array([
            [m[0], m[1], m[2], m[3]],
            [m[1], m[2], m[3], m[4]],
            [m[2], m[3], m[4], m[5]],
            [m[3], m[4], m[5], m[6]]
        ])
        
        # M_ij = <x^{i+j}>, i,j from {1,3,5,7}
        M_odd = np.array([
            [m[1], m[2], m[3], m[4]],
            [m[2], m[3], m[4], m[5]],
            [m[3], m[4], m[5], m[6]],
            [m[4], m[5], m[6], m[7]]
        ])
        
        # A matrix is positive semidefinite if all its eigenvalues are non-negative.
        # Use a small tolerance for floating point inaccuracies.
        tol = -1e-9 
        try:
            # eigvalsh is for symmetric/Hermitian matrices and is faster
            eig_even = np.linalg.eigvalsh(M_even)
            eig_odd = np.linalg.eigvalsh(M_odd)
        except np.linalg.LinAlgError:
            return False
        
        return np.all(eig_even >= tol) and np.all(eig_odd >= tol)

    # Grid search to find the minimal E and corresponding m1
    # Literature values suggest E is around 1.3 and m1 around 0.4
    E_values = np.arange(1.3, 1.4, 0.001)
    m1_values = np.arange(0.4, 0.5, 0.001)
    
    min_E_found = float('inf')
    min_m1_at_min_E = float('inf')
    
    for E in E_values:
        found_valid_m1_for_this_E = False
        current_min_m1_for_this_E = float('inf')
        
        for m1 in m1_values:
            if check_psd(E, m1):
                found_valid_m1_for_this_E = True
                # Since we iterate m1 from low to high, the first one we find is the minimum
                current_min_m1_for_this_E = m1
                break  # Found the smallest valid m1 for this E, move to next E
        
        if found_valid_m1_for_this_E:
            # Since we iterate E from low to high, the first E that works is the minimum
            min_E_found = E
            min_m1_at_min_E = current_min_m1_for_this_E
            break # Found the ground state, exit the search

    if min_E_found == float('inf'):
        print("No solution found in the specified range.")
        print("You may need to expand the search grid for E and/or <x^2>.")
    else:
        # We need to output each number in the final result.
        e_val_str = f"{min_E_found:.3f}"
        x2_val_str = f"{min_m1_at_min_E:.3f}"
        
        print(f"For V(x) = x^2 + x^4 and K = 7:")
        print(f"Minimal Energy E = {e_val_str}")
        print(f"Corresponding <x^2> = {x2_val_str}")
        print("<<<Minimal E is {} and minimal <x^2> is {}>>>".format(e_val_str, x2_val_str))


solve_bootstrap_qm()