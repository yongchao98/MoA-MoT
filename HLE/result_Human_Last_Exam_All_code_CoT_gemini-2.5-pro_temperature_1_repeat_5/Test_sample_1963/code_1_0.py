import numpy as np

def calculate_moments(E, x2_val):
    """
    Calculates the first 8 even moments m_k = <x^{2k}> for k=0 to 7.
    The calculation is based on the recursion relation derived for V(x) = x^2 + x^4.
    
    Args:
        E (float): The energy parameter.
        x2_val (float): The value of <x^2> (which is m_1).

    Returns:
        numpy.ndarray or None: An array of moments m_0 to m_7, or None if a physical
                               constraint (like m_k >= 0) is obviously violated early.
    """
    m = np.zeros(8)
    m[0] = 1.0
    m[1] = x2_val

    # The recursion is derived from the prompt for t=2k+1.
    # (8k+12)m_{k+2} = (2k+1)(2k)(2k-1)m_{k-1} + 4(2k+1)E m_k - (8k+8)m_{k+1}
    
    # For k=0 (t=1): calculates m_2 = <x^4>
    # 12 * m_2 = 4*E*m_0 - 8*m_1
    m[2] = (4 * E * m[0] - 8 * m[1]) / 12.0
    # A necessary condition is <x^4> >= 0.
    if m[2] < 0:
        return None

    # For k=1 (t=3): calculates m_3 = <x^6>
    # 20 * m_3 = 6*m_0 + 12*E*m_1 - 16*m_2
    m[3] = (6 * m[0] + 12 * E * m[1] - 16 * m[2]) / 20.0
    
    # For k=2 (t=5): calculates m_4 = <x^8>
    # 28 * m_4 = 60*m_1 + 20*E*m_2 - 24*m_3
    m[4] = (60 * m[1] + 20 * E * m[2] - 24 * m[3]) / 28.0

    # For k=3 (t=7): calculates m_5 = <x^10>
    # 36 * m_5 = 210*m_2 + 28*E*m_3 - 32*m_4
    m[5] = (210 * m[2] + 28 * E * m[3] - 32 * m[4]) / 36.0

    # For k=4 (t=9): calculates m_6 = <x^12>
    # 44 * m_6 = 504*m_3 + 36*E*m_4 - 40*m_5
    m[6] = (504 * m[3] + 36 * E * m[4] - 40 * m[5]) / 44.0

    # For k=5 (t=11): calculates m_7 = <x^14>
    # 52 * m_7 = 990*m_4 + 44*E*m_5 - 48*m_6
    m[7] = (990 * m[4] + 44 * E * m[5] - 48 * m[6]) / 52.0
    
    return m

def check_positive_semidefinite(moments):
    """
    Checks if the two moment matrices M_e and M_o are positive semidefinite.
    A matrix is positive semidefinite if all its eigenvalues are non-negative.
    
    Args:
        moments (numpy.ndarray): The array of even moments m_0 to m_7.

    Returns:
        bool: True if both matrices are positive semidefinite, False otherwise.
    """
    m = moments
    
    # M_e is constructed from moments <x^{2i+2j}> = m_{i+j} for i,j in {0,1,2,3}
    M_e = np.array([
        [m[0], m[1], m[2], m[3]],
        [m[1], m[2], m[3], m[4]],
        [m[2], m[3], m[4], m[5]],
        [m[3], m[4], m[5], m[6]]
    ])
    
    # M_o is constructed from moments <x^{(2i+1)+(2j+1)}> = m_{i+j+1} for i,j in {0,1,2,3}
    M_o = np.array([
        [m[1], m[2], m[3], m[4]],
        [m[2], m[3], m[4], m[5]],
        [m[3], m[4], m[5], m[6]],
        [m[4], m[5], m[6], m[7]]
    ])
    
    # Use a small tolerance for floating point inaccuracies
    tolerance = -1e-9
    
    # np.linalg.eigvalsh is used for symmetric/Hermitian matrices and is faster.
    eigenvalues_e = np.linalg.eigvalsh(M_e)
    if np.min(eigenvalues_e) < tolerance:
        return False
        
    eigenvalues_o = np.linalg.eigvalsh(M_o)
    if np.min(eigenvalues_o) < tolerance:
        return False
        
    return True

def find_minimal_values():
    """
    Performs a grid search to find the minimal E and corresponding <x^2>.
    """
    # Define search ranges based on known approximations for this potential.
    # This makes the search more efficient.
    # The search step of 0.001 is chosen to achieve the desired 3-digit precision.
    e_range = np.arange(1.3, 1.45, 0.001)
    x2_range = np.arange(0.3, 0.45, 0.001)

    min_E_found = float('inf')
    best_x2_found = -1

    # Iterate E from low to high. The first E that allows a valid <x^2> is the minimum.
    for E in e_range:
        found_valid_solution = False
        for x2_val in x2_range:
            moments = calculate_moments(E, x2_val)
            
            # If moments are invalid (e.g., negative <x^4>), skip.
            if moments is None:
                continue
            
            # Check the positivity constraint.
            if check_positive_semidefinite(moments):
                min_E_found = E
                best_x2_found = x2_val
                found_valid_solution = True
                # We found the first valid pair for this E, so we can break the inner loop.
                # And since we iterate E upwards, this is the minimal E.
                break
        
        if found_valid_solution:
            break
            
    # Output the final results
    print("The minimal value of E and the corresponding <x^2> are:")
    print(f"Minimal E = {min_E_found:.3f}")
    print(f"Corresponding <x^2> = {best_x2_found:.3f}")

if __name__ == '__main__':
    find_minimal_values()