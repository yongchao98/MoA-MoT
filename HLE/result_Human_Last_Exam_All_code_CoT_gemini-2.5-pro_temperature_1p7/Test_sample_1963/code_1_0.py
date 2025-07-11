import numpy as np

def calculate_moments(E, m1, max_k_index):
    """
    Calculates the even moments m_k = <x^{2k}> using the recursion relation
    for the potential V(x) = x^2 + x^4.

    Args:
        E (float): The energy value to test.
        m1 (float): The value for <x^2> to test.
        max_k_index (int): The maximum index k for the moment m_k to calculate.

    Returns:
        numpy.ndarray: An array containing moments from m_0 to m_{max_k_index}.
    """
    m = np.zeros(max_k_index + 1)
    if max_k_index >= 0:
        m[0] = 1.0  # <x^0> = 1
    if max_k_index >= 1:
        m[1] = m1   # <x^2> is a test parameter

    # The recursion relation is:
    # 4(2k+1) m_{k+1} = 4(2k-1)E m_{k-1} + (2k-1)(2k-2)(2k-3) m_{k-2} - 8k m_k
    # This relation is valid for k >= 1. We calculate m[k+1] in the loop.
    for k in range(1, max_k_index):
        # Retrieve lower-order moments from the array. m_k corresponds to m[k].
        m_k = m[k]
        m_k_minus_1 = m[k-1]
        
        # For k=1, the m_{k-2} term (m_{-1}) is zero as per the bootstrap rules.
        m_k_minus_2 = m[k-2] if k > 1 else 0.0

        # Apply the recursion formula
        term1 = 4 * (2 * k - 1) * E * m_k_minus_1
        term2 = (2 * k - 1) * (2 * k - 2) * (2 * k - 3) * m_k_minus_2
        term3 = -8 * k * m_k
        
        numerator = term1 + term2 + term3
        denominator = 4 * (2 * k + 1)
        
        # Handle potential division by zero, though not expected for k>=1
        if denominator == 0:
            print("Error: Division by zero in recursion.")
            return None
        
        m[k + 1] = numerator / denominator
        
    return m

def is_positive_semidefinite(matrix):
    """
    Checks if a matrix is positive semidefinite by verifying that all its
    eigenvalues are non-negative. It assumes the matrix is symmetric.
    """
    # Use eigvalsh for symmetric matrices for efficiency and numerical stability
    eigenvalues = np.linalg.eigvalsh(matrix)
    # Allow for small numerical errors
    return np.all(eigenvalues >= -1e-9)

def find_minimal_values():
    """
    Performs the grid search to find the minimal E and <x^2>.
    """
    # K=7 requires an 8x8 moment matrix M_{ij} = <x^{i+j}>
    # This matrix decouples into two 4x4 matrices.
    # M_even requires moments m_0...m_6.
    # M_odd requires moments m_1...m_7.
    # So, we need to calculate moments up to m_7.
    max_k_index = 7

    # Define a search grid for E and <x^2>.
    # Based on known results, the ground state energy is around 1.39.
    e_range = np.arange(1.3, 1.5, 0.001)
    x2_range = np.arange(0.4, 0.7, 0.001)

    min_E = float('inf')
    best_x2 = -1
    
    # Iterate through the grid, starting with the lowest E
    for e_val in e_range:
        for x2_val in x2_range:
            # 1. Calculate moments for the current (E, <x^2>) pair
            moments = calculate_moments(e_val, x2_val, max_k_index)

            # 2. Construct the positivity matrices
            # M_even is a 4x4 Hankel matrix of moments m_0 to m_6
            M_even = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    M_even[i, j] = moments[i + j]

            # M_odd is a 4x4 Hankel matrix of moments m_1 to m_7
            M_odd = np.zeros((4, 4))
            for i in range(4):
                for j in range(4):
                    M_odd[i, j] = moments[i + j + 1]

            # 3. Check if both matrices are positive semidefinite
            if is_positive_semidefinite(M_even) and is_positive_semidefinite(M_odd):
                min_E = e_val
                best_x2 = x2_val
                
                # The first valid pair found gives the minimal E and corresponding minimal <x^2>
                # because we are iterating E and <x^2> in increasing order.
                print("Found a valid solution at the boundary of the allowed region.")
                print("The final values are determined by this boundary condition.")
                print(f"Minimal E: {min_E:.3f}")
                print(f"Corresponding minimal <x^2>: {best_x2:.3f}")
                return min_E, best_x2

    print("No solution found in the specified search range.")
    return None, None

if __name__ == "__main__":
    E, x2 = find_minimal_values()
    if E is not None:
        print("\n--- Final Answer ---")
        # The prompt expects a single numerical value in the final tag.
        # Typically, the ground state energy E is the primary quantity of interest.
        print(f"The minimal ground state energy E is {E:.3f}")
        print(f"<<<{E:.3f}>>>")
