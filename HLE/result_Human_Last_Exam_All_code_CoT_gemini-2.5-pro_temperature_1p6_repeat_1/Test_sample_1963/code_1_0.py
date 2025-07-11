import numpy as np

def calculate_moments(E, m1):
    """
    Calculates the even moments m_k = <x^{2k}> using the recursion relations.
    We need moments m_0 through m_7.
    """
    m = np.zeros(8)
    
    # Initial conditions
    m[0] = 1.0  # <x^0>
    m[1] = m1   # <x^2>
    
    # Calculate m_2 from the t=1 relation: <x^4> = E/3 - (2/3)<x^2>
    m[2] = (E / 3.0) - (2.0 / 3.0) * m[1]
    
    # Use the general three-term recurrence relation for k >= 1:
    # (8k+12)m_{k+2} = (2k+1)(2k)(2k-1)m_{k-1} + 4(2k+1)Em_k - (8k+8)m_{k+1}
    for k in range(1, 6):  # Loop to calculate m_3, m_4, m_5, m_6, m_7
        numerator = (
            (2*k + 1) * (2*k) * (2*k - 1) * m[k - 1]
            + 4 * (2*k + 1) * E * m[k]
            - (8*k + 8) * m[k + 1]
        )
        denominator = 8*k + 12
        if denominator == 0:
            return None # Avoid division by zero
        m[k + 2] = numerator / denominator
        
    return m

def is_psd(matrix, tol=1e-9):
    """Checks if a matrix is positive semidefinite by checking its eigenvalues."""
    # Use eigvalsh for symmetric matrices; it's faster and more stable.
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.all(eigenvalues >= -tol)

def is_feasible(E, m1, tol=1e-9):
    """Checks if a given (E, <x^2>) pair is physically allowed."""
    
    # Fundamental constraint from Cauchy-Schwarz: <x^4> >= <x^2>^2, so m_2 >= m_1^2
    m2_val = (E / 3.0) - (2.0 / 3.0) * m1
    if m2_val < m1**2 - tol:
        return False
        
    moments = calculate_moments(E, m1)
    if moments is None or np.any(np.isnan(moments)):
        return False
        
    # All even moments must be non-negative
    if np.any(moments < -tol):
        return False

    # Construct the Hankel matrices M_E and M_O for K=7
    # M_E uses moments m_0 to m_6
    # M_O uses moments m_1 to m_7
    M_E = np.zeros((4, 4))
    M_O = np.zeros((4, 4))
    
    for i in range(4):
        for j in range(4):
            M_E[i, j] = moments[i + j]
            M_O[i, j] = moments[i + j + 1]

    # Check if both matrices are positive semidefinite
    return is_psd(M_E, tol) and is_psd(M_O, tol)

def solve():
    """
    Performs a grid search to find the minimal E and <x^2>.
    """
    # Grid search parameters to ensure 3-digit precision
    e_grid = np.arange(1.37, 1.39, 0.001)
    x2_grid = np.arange(0.47, 0.49, 0.001)

    min_E_found = float('inf')
    
    # Find the minimum allowed energy E
    for e_val in e_grid:
        for x2_val in x2_grid:
            if is_feasible(e_val, x2_val):
                min_E_found = e_val
                break # Found the first (minimal) E, stop searching E
        if min_E_found != float('inf'):
            break
            
    if min_E_found == float('inf'):
        print("No solution found in the specified grid. Please expand the search range.")
        return

    # For the minimal E, find the minimal corresponding <x^2>
    min_x2_found = float('inf')
    for x2_val in x2_grid:
        if is_feasible(min_E_found, x2_val):
            min_x2_found = x2_val
            break # Found the first (minimal) x2, stop
    
    print(f"Minimal E = {min_E_found:.3f}")
    print(f"Minimal <x^2> = {min_x2_found:.3f}\n")
    
    # Display the final moment matrices for these values
    final_moments = calculate_moments(min_E_found, min_x2_found)
    M_E = np.zeros((4, 4))
    M_O = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            M_E[i, j] = final_moments[i + j]
            M_O[i, j] = final_moments[i + j + 1]
    
    np.set_printoptions(precision=5, suppress=True)
    print("Resulting Moment Matrix for Even operators (M_E):")
    print(M_E)
    print("\nResulting Moment Matrix for Odd operators (M_O):")
    print(M_O)
    
    # Print final answer in the requested format
    print(f"\n<<<E={min_E_found:.3f}, <x^2>={min_x2_found:.3f}>>>")

# Run the solver
solve()