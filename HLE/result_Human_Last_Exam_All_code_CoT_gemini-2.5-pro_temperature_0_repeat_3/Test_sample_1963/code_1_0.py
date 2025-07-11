import numpy as np

def calculate_moments(E, x2_val, max_n=7):
    """
    Calculates the even moments m_k = <x^{2k}> using the recursion relation.
    The recursion is derived for V(x) = x^2 + x^4.
    We need moments up to m_7 for K=7.
    """
    m = np.zeros(max_n + 1)
    
    # Initial conditions
    m[0] = 1.0  # <x^0>
    if max_n >= 1:
        m[1] = x2_val  # <x^2>

    # The first relation from t=1 gives m_2 = <x^4>
    # 4E - 8<x^2> - 12<x^4> = 0
    if max_n >= 2:
        m[2] = (4 * E * m[0] - 8 * m[1]) / 12.0

    # General recursion for m_{k+2} from k=1 onwards
    # (8k+12)m_{k+2} = (2k+1)(2k)(2k-1)m_{k-1} + 4(2k+1)E m_k - (8k+8)m_{k+1}
    for k in range(1, max_n - 1):
        if k + 2 > max_n:
            break
        
        numerator = ((2*k + 1) * (2*k) * (2*k - 1) * m[k-1] +
                     4 * (2*k + 1) * E * m[k] -
                     (8*k + 8) * m[k+1])
        denominator = 8*k + 12
        
        m[k+2] = numerator / denominator
        
    return m

def is_positive_semidefinite(matrix):
    """
    Checks if a symmetric matrix is positive semidefinite by checking its eigenvalues.
    A small tolerance is used for floating point inaccuracies.
    """
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.all(eigenvalues >= -1e-9)

def check_constraints(E, x2_val):
    """
    For a given E and <x^2>, calculates the moments and checks if the
    moment matrices M_e and M_o are positive semidefinite.
    """
    if x2_val < 0:
        return False
        
    # For K=7, we need moments up to m_7 = <x^14>
    moments = calculate_moments(E, x2_val, max_n=7)
    m = moments
    
    # Construct the even-power matrix M_e
    # M_e_{ij} = m_{i+j} for i,j in {0,1,2,3}
    M_e = np.array([
        [m[0], m[1], m[2], m[3]],
        [m[1], m[2], m[3], m[4]],
        [m[2], m[3], m[4], m[5]],
        [m[3], m[4], m[5], m[6]]
    ])
    
    # Construct the odd-power matrix M_o
    # M_o_{ij} = m_{i+j+1} for i,j in {0,1,2,3}
    M_o = np.array([
        [m[1], m[2], m[3], m[4]],
        [m[2], m[3], m[4], m[5]],
        [m[3], m[4], m[5], m[6]],
        [m[4], m[5], m[6], m[7]]
    ])

    # Check if both matrices are positive semidefinite
    return is_positive_semidefinite(M_e) and is_positive_semidefinite(M_o)

def find_minimal_values():
    """
    Performs a numerical search to find the minimal E and corresponding <x^2>.
    """
    # Step 1: Binary search for the minimal allowed energy E.
    E_low = 1.0
    E_high = 2.0
    min_E_candidate = float('inf')

    # 25 iterations for high precision
    for _ in range(25):
        E_mid = (E_low + E_high) / 2
        is_allowed = False
        # Check if any x2 value is allowed for this E
        x2_search_range = np.linspace(0.1, 0.8, 151)
        for x2 in x2_search_range:
            if check_constraints(E_mid, x2):
                is_allowed = True
                break
        
        if is_allowed:
            # This E is allowed, so it's a candidate for the minimum.
            # The true minimum could be even lower.
            min_E_candidate = E_mid
            E_high = E_mid
        else:
            # This E is too low, no physical state exists.
            E_low = E_mid

    final_E = min_E_candidate

    # Step 2: For the found minimal E, find the corresponding <x^2>.
    # At the minimal E, the allowed region for <x^2> should shrink to a point.
    allowed_x2_vals = []
    # Use a very fine grid for x2 to find the narrow allowed range
    x2_search_range_final = np.linspace(0.1, 0.8, 1001)
    for x2 in x2_search_range_final:
        # We check at a slightly higher E to ensure the interval is not empty due to precision limits
        if check_constraints(final_E * 1.000001, x2):
            allowed_x2_vals.append(x2)

    if not allowed_x2_vals:
        # Should not happen with the slight E increase, but as a fallback
        return -1, -1

    final_x2 = np.mean(allowed_x2_vals)
    
    return final_E, final_x2

if __name__ == '__main__':
    min_E, min_x2 = find_minimal_values()
    
    print("Based on the bootstrap method with K=7:")
    print(f"The minimal value of E is: {min_E:.3f}")
    print(f"The corresponding value of <x^2> is: {min_x2:.3f}")
    
    # Final answer in the required format
    print(f"<<<{min_x2:.3f}, {min_E:.3f}>>>")
