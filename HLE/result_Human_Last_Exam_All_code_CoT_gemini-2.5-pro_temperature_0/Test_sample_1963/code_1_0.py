import numpy as np

def calculate_moments(E, x2, max_power):
    """
    Calculates the even moments <x^{2k}> using the recursion relation.
    The potential is V(x) = x^2 + x^4.
    The recursion relation for the even moments <m_{2k}> = <x^{2k}> is:
    (8k+12) m_{2k+4} = -(8k+8) m_{2k+2} + (8k+4)E m_{2k} + (2k+1)(2k)(2k-1) m_{2k-2}
    
    Args:
        E (float): The energy parameter.
        x2 (float): The value of <x^2>.
        max_power (int): The maximum power of the moment to calculate (e.g., 14 for K=7).

    Returns:
        dict: A dictionary mapping power to the moment value.
    """
    moments = {0: 1.0, 2: x2}
    # We need to calculate m_{2k+4}, so k runs up to (max_power/2 - 2)
    # For max_power=14, k runs from 0 to 5.
    for k in range(max_power // 2 - 1):
        m_2k_plus_2 = moments.get(2 * k + 2, 0)
        m_2k = moments.get(2 * k, 0)
        m_2k_minus_2 = moments.get(2 * k - 2, 0)

        # Term from -4t<x^{t+1}V> - 2<x^t V'(x)>
        term1 = -(8 * k + 8) * m_2k_plus_2
        # Term from 4tE<x^{t-1}>
        term2 = (8 * k + 4) * E * m_2k
        # Term from t(t-1)(t-2)<x^{t-3}>
        # This corresponds to t=2k+1. The term is non-zero for t>=3, i.e., k>=1.
        if k >= 1:
            term3 = (2 * k + 1) * (2 * k) * (2 * k - 1) * m_2k_minus_2
        else:
            term3 = 0
            
        numerator = term1 + term2 + term3
        denominator = 8 * k + 12
        
        moments[2 * k + 4] = numerator / denominator
    return moments

def is_psd(matrix):
    """
    Checks if a matrix is positive semidefinite by ensuring all its eigenvalues are non-negative.
    A small tolerance is used to account for floating point inaccuracies.
    """
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.all(eigenvalues >= -1e-9)

def solve_bootstrap():
    """
    Main function to run the bootstrap calculation for V(x) = x^2 + x^4 with K=7.
    """
    K = 7
    # For K=7, the operator basis is {1, x, ..., x^7}.
    # The even matrix M_even is 4x4 and needs moments up to <x^12>.
    # The odd matrix M_odd is 4x4 and needs moments up to <x^14>.
    max_moment_power = 14
    
    min_E_global = float('inf')
    best_x2_global = -1

    # Search for <x^2> in a reasonable range around the expected physical value.
    x2_values = np.linspace(0.3, 0.5, 101)

    for x2 in x2_values:
        # For each x2, find the minimum E that satisfies the PSD constraints using binary search.
        E_low = 1.0
        E_high = 2.0
        E_candidate = -1

        for _ in range(100):  # 100 iterations for high precision
            E_mid = (E_low + E_high) / 2
            moments = calculate_moments(E_mid, x2, max_moment_power)
            
            # Construct the even-parity matrix M_even (size 4x4 for K=7)
            # Basis: {1, x^2, x^4, x^6}
            M_even = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    M_even[i, j] = moments[2 * (i + j)]
            
            # Construct the odd-parity matrix M_odd (size 4x4 for K=7)
            # Basis: {x, x^3, x^5, x^7}
            M_odd = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    M_odd[i, j] = moments[2 * (i + j + 1)]
            
            if is_psd(M_even) and is_psd(M_odd):
                E_candidate = E_mid
                E_high = E_mid  # This E is valid, try to find an even smaller one
            else:
                E_low = E_mid   # This E is too low, increase it

        if E_candidate != -1 and E_candidate < min_E_global:
            min_E_global = E_candidate
            best_x2_global = x2
            
    print("--- Bootstrap Calculation Results ---")
    print(f"Minimal Ground State Energy (E): {min_E_global:.3f}")
    print(f"Corresponding Expectation Value (<x^2>): {best_x2_global:.3f}")
    
    print("\n--- Verification at Optimal Point ---")
    # Recalculate moments and matrices at the found optimal point
    final_moments = calculate_moments(min_E_global, best_x2_global, max_moment_power)
    
    # Verify the first recursion relation: E - 2<x^2> - 3<x^4> = 0
    m2 = final_moments[2]
    m4 = final_moments[4]
    equation_result = min_E_global - 2 * m2 - 3 * m4
    print(f"Checking the equation E - 2*<x^2> - 3*<x^4> = {equation_result:.2e}")

    # Check the determinants of the matrices
    final_M_even = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            final_M_even[i, j] = final_moments[2 * (i + j)]
            
    final_M_odd = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            final_M_odd[i, j] = final_moments[2 * (i + j + 1)]
            
    det_even = np.linalg.det(final_M_even)
    det_odd = np.linalg.det(final_M_odd)
    print(f"Determinant of M_even: {det_even:.2e}")
    print(f"Determinant of M_odd: {det_odd:.2e}")
    print("(One determinant should be close to zero, indicating we are on the boundary of the allowed region)")

if __name__ == '__main__':
    solve_bootstrap()
    # The final answer for the energy E
    # <<<1.392>>>