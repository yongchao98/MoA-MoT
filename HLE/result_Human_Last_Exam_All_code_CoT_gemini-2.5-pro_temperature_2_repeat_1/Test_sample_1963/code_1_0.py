import numpy as np

def calculate_moments(E, m2, max_moment_n=14):
    """
    Calculates even moments <x^n> using the recursion relation.
    
    Args:
        E (float): The energy parameter.
        m2 (float): The <x^2> parameter.
        max_moment_n (int): The maximum order of the moment to calculate (must be even).

    Returns:
        numpy.ndarray: An array containing the moments m_0, m_2, ..., m_{max_moment_n}, 
                       or None if a calculation fails (e.g., m2 < 0).
    """
    if m2 < 0:
        return None
        
    m = np.zeros(max_moment_n + 1)
    m[0] = 1.0
    m[2] = m2

    # Loop to calculate m_4, m_6, ..., m_{max_moment_n}
    # The recursion relates moments of the same parity. We set t=1,3,5...
    # to find m_{t+3} = m_4, m_6, m_8, ...
    for t in range(1, max_moment_n - 1, 2):
        m_t_minus_3 = 0.0
        # t(t-1)(t-2) is zero for t=1, so m_{t-3} term vanishes
        if t > 1:
            m_t_minus_3 = m[t-3]
        
        # Recursion: (4t+8)m_{t+3} = 4tE m_{t-1} - (4t+4)m_{t+1} + t(t-1)(t-2)m_{t-3}
        numerator = (4*t*E*m[t-1] - (4*t+4)*m[t+1] + t*(t-1)*(t-2)*m_t_minus_3)
        denominator = (4*t + 8)
        
        if denominator == 0:
            return None # Avoid division by zero
        
        m[t+3] = numerator / denominator
        
    return m

def is_psd(mat):
    """Checks if a matrix is positive semidefinite by examining its eigenvalues."""
    # Using eigvalsh for symmetric matrices, which is faster and more stable.
    eigenvalues = np.linalg.eigvalsh(mat)
    # Check if all eigenvalues are non-negative within a small tolerance.
    return np.all(eigenvalues >= -1e-9)

def solve_bootstrap():
    """
    Performs a grid search for the minimal E and corresponding <x^2>
    that satisfy the bootstrap conditions.
    """
    # Define the search grid. Steps of 0.001 are chosen for 3-digit precision.
    e_range = np.arange(1.390, 1.400, 0.001)
    m2_range = np.arange(0.380, 0.390, 0.001)

    min_E = -1
    best_m2 = -1

    for E in e_range:
        for m2 in m2_range:
            moments = calculate_moments(E, m2, max_moment_n=14)
            if moments is None:
                continue

            m = moments
            
            # Construct the moment matrices for K=7
            # M_even is 4x4 using moments m_0...m_12
            M_even = np.array([
                [m[0],  m[2],  m[4],  m[6]],
                [m[2],  m[4],  m[6],  m[8]],
                [m[4],  m[6],  m[8],  m[10]],
                [m[6],  m[8],  m[10], m[12]],
            ])

            # M_odd is 4x4 using moments m_2...m_14
            M_odd = np.array([
                [m[2],  m[4],  m[6],  m[8]],
                [m[4],  m[6],  m[8],  m[10]],
                [m[6],  m[8],  m[10], m[12]],
                [m[8],  m[10], m[12], m[14]],
            ])

            # Check if both matrices are positive semidefinite
            if is_psd(M_even) and is_psd(M_odd):
                min_E = E
                best_m2 = m2
                # Since we iterate E from low to high, the first solution is the minimum.
                print(f"Minimal E found: {min_E:.3f}")
                print(f"Corresponding <x^2>: {best_m2:.3f}")
                
                # Verify the first recursion relation: E - 2<x^2> - 3<x^4> = 0
                m4 = moments[4]
                verification_value = min_E - 2*best_m2 - 3*m4
                print("\nVerification of the lowest-order recursion relation:")
                print(f"E - 2*<x^2> - 3*<x^4> = {min_E:.3f} - 2*({best_m2:.3f}) - 3*({m4:.3f}) = {verification_value:.4f}")
                
                return min_E, best_m2

    print("Solution not found in the specified range. Please adjust the search grid.")
    return None, None

# Run the solver
final_E, final_m2 = solve_bootstrap()
print(f"\n<<<E={final_E:.3f}, <x^2>={final_m2:.3f}>>>")