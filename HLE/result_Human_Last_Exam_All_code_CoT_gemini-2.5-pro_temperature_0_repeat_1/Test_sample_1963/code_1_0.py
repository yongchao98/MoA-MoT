import numpy as np

def calculate_moments(E, m1, K):
    """
    Calculates the moments <x^{2k}> for k=0 to K using the recursion relation.
    V(x) = x^2 + x^4
    """
    # For K=7, we need moments m_0 to m_7.
    num_moments = K + 1
    moments = [0.0] * num_moments

    # m_0 = <x^0> = 1
    moments[0] = 1.0
    # m_1 = <x^2> is given
    moments[1] = m1

    # From <p^2> + <x^2> + <x^4> = E and <p^2> = <x^2> + 2<x^4>, we get 2<x^2> + 3<x^4> = E.
    # So, m_2 = <x^4> = (E - 2*m_1) / 3
    m2 = (E - 2.0 * m1) / 3.0
    # <x^4> must be non-negative.
    if m2 < 0:
        return None
    moments[2] = m2

    # The recursion relation for m_k = <x^{2k}> is:
    # m_{k+2} = (1/(8k+12)) * [ (2k+1)(2k)(2k-1)m_{k-1} + 4(2k+1)E*m_k - (8k+8)m_{k+1} ]
    # We need to calculate up to m_7, so we loop for k from 1 to 5.
    for k in range(1, K - 1): # k=1..5, this calculates m_3..m_7
        m_k_minus_1 = moments[k-1]
        m_k = moments[k]
        m_k_plus_1 = moments[k+1]

        term1 = (2*k+1)*(2*k)*(2*k-1) * m_k_minus_1
        term2 = 4*(2*k+1)*E * m_k
        term3 = -(8*k+8) * m_k_plus_1
        
        numerator = term1 + term2 + term3
        denominator = 8*k + 12
        
        if denominator == 0: return None # Should not happen for k>=1
        
        m_k_plus_2 = numerator / denominator
        moments[k+2] = m_k_plus_2
        
    return moments

def check_psd(moments, K):
    """
    Constructs the M_even and M_odd matrices and checks if they are positive semidefinite.
    """
    # For K=7, both matrices are 4x4.
    # M_even_{ij} = m_{i+j} for i,j in {0,1,2,3}
    size_even = K // 2 + 1
    M_even = np.zeros((size_even, size_even))
    for i in range(size_even):
        for j in range(size_even):
            M_even[i, j] = moments[i+j]
            
    # M_odd_{ij} = m_{i+j+1} for i,j in {0,1,2,3}
    size_odd = (K - 1) // 2 + 1
    M_odd = np.zeros((size_odd, size_odd))
    for i in range(size_odd):
        for j in range(size_odd):
            M_odd[i, j] = moments[i+j+1]

    # Check positive semidefiniteness using eigenvalues.
    # A small negative tolerance is used for numerical stability.
    eig_even = np.linalg.eigvalsh(M_even)
    if not np.all(eig_even >= -1e-9):
        return False
        
    eig_odd = np.linalg.eigvalsh(M_odd)
    if not np.all(eig_odd >= -1e-9):
        return False
        
    return True

def solve_bootstrap():
    """
    Performs a grid search to find the minimal E and corresponding <x^2>.
    """
    K = 7
    
    # Based on known results, the ground state energy is around 1.39.
    # We search in a narrow grid to find the solution with 3-digit precision.
    E_range = np.arange(1.390, 1.400, 0.001)
    x2_range = np.arange(0.380, 0.400, 0.001)
    
    found_E = None
    found_x2 = None

    # Iterate through the grid. The first valid (E, <x^2>) pair will be the minimum.
    for E in E_range:
        for m1 in x2_range:
            if m1 <= 0:
                continue
            
            moments = calculate_moments(E, m1, K)
            
            if moments is None:
                continue
            
            if check_psd(moments, K):
                found_E = E
                found_x2 = m1
                # Break both loops once the first valid pair is found
                break
        if found_E is not None:
            break

    if found_E is None:
        print("No solution found. Please consider expanding the search ranges for E and <x^2>.")
        return

    print("Found minimal ground state energy and expectation value:")
    print(f"E = {found_E:.3f}")
    print(f"<x^2> = {found_x2:.3f}")
    print("-" * 40)

    # For the final answer, show the numbers in the final constraint equations.
    final_moments = calculate_moments(found_E, found_x2, K)
    
    print("The final equations are the positive semidefinite constraints on the moment matrices.")
    print("Calculated moments <x^{2n}> for n=0 to 7:")
    for i, m in enumerate(final_moments):
        print(f"m_{i} = <x^{2*i}> = {m:.4f}")
    print("-" * 40)

    # Construct and print final matrices and their eigenvalues to verify the solution
    size_even = K // 2 + 1
    M_even = np.array([[final_moments[i+j] for j in range(size_even)] for i in range(size_even)])

    size_odd = (K - 1) // 2 + 1
    M_odd = np.array([[final_moments[i+j+1] for j in range(size_odd)] for i in range(size_odd)])
            
    print("Final matrix M_even (from <x^{2i}x^{2j}>):")
    print(np.round(M_even, 4))
    eig_even = np.linalg.eigvalsh(M_even)
    print("\nEigenvalues of M_even (all non-negative):")
    print(np.round(eig_even, 6))
    print("-" * 40)

    print("Final matrix M_odd (from <x^{2i+1}x^{2j+1}>):")
    print(np.round(M_odd, 4))
    eig_odd = np.linalg.eigvalsh(M_odd)
    print("\nEigenvalues of M_odd (all non-negative):")
    print(np.round(eig_odd, 6))
    print("-" * 40)

# Run the solver
solve_bootstrap()