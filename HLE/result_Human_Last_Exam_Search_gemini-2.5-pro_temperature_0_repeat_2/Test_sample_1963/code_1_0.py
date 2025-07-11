import numpy as np

def calculate_moments(E, A, K):
    """Calculates even moments <x^n> up to n=2K+2."""
    max_n = 2 * (K + 1) # For K=7, we need up to x^14 for M_odd
    moments = {0: 1.0, 2: A}
    
    # The potential is symmetric, so odd moments are zero.
    # We only need to calculate even moments.
    # The recursion relation connects moments with powers differing by 2.
    # Let's use the version from the thought process:
    # <x^{k+4}> = (4(k+1)E <x^k> - (4k+8) <x^{k+2}> + (k+1)k(k-1) <x^{k-2}>) / (4k+12)
    
    for k in range(0, max_n, 2):
        if k + 4 > max_n:
            break
        
        term3 = 0
        if k >= 2:
            term3 = (k + 1) * k * (k - 1) * moments[k - 2]
            
        numerator = 4 * (k + 1) * E * moments[k] - (4 * k + 8) * moments[k + 2] + term3
        denominator = 4 * k + 12
        
        if denominator == 0:
            # This should not happen for k>=0
            return None
            
        moments[k + 4] = numerator / denominator
        
    return moments

def check_psd(E, A, K):
    """Checks if the moment matrices are positive semidefinite."""
    if A < 0:
        return False
        
    moments = calculate_moments(E, A, K)
    if moments is None:
        return False

    # For K=7, the operators are {1, x, x^2, ..., x^7}.
    # The even-power operator is O_even = c0 + c2*x^2 + c4*x^4 + c6*x^6.
    # The odd-power operator is O_odd = c1*x + c3*x^3 + c5*x^5 + c7*x^7.
    # This leads to two 4x4 matrices.
    
    num_even = (K // 2) + 1 if K % 2 == 1 else K // 2 + 1 # 4 for K=7
    num_odd = (K // 2) + 1 if K % 2 == 1 else K // 2 # 4 for K=7

    # M_even: <x^{2i+2j}> for i,j in 0..3
    M_even = np.zeros((num_even, num_even))
    for i in range(num_even):
        for j in range(num_even):
            power = 2 * (i + j)
            if power not in moments: return False
            M_even[i, j] = moments[power]

    # M_odd: <x^{2i+2j+2}> for i,j in 0..3
    M_odd = np.zeros((num_odd, num_odd))
    for i in range(num_odd):
        for j in range(num_odd):
            power = 2 * (i + j) + 2
            if power not in moments: return False
            M_odd[i, j] = moments[power]

    # Check if matrices are positive semidefinite by checking eigenvalues
    # A small tolerance is used to account for floating point inaccuracies.
    if not np.all(np.linalg.eigvalsh(M_even) >= -1e-9):
        return False
    if not np.all(np.linalg.eigvalsh(M_odd) >= -1e-9):
        return False
        
    return True

# Grid search to find the minimal E
def find_minimums(K):
    # Search range based on literature and preliminary analysis
    E_range = np.linspace(1.38, 1.40, 101)
    A_range = np.linspace(0.4, 0.5, 101)

    min_E_found = float('inf')
    corresponding_A = -1

    for E in E_range:
        for A in A_range:
            if check_psd(E, A, K):
                # Since we are iterating E from low to high, the first E
                # that has a valid A is the minimum E.
                min_E_found = E
                corresponding_A = A
                # We can stop searching for this E and move to the next one
                # to map the boundary, but for finding the minimum, we can stop here.
                # To get the most precise value from the grid, we find the first valid (E,A) pair.
                return min_E_found, corresponding_A
    
    return min_E_found, corresponding_A

# Set K=7 and run the search
K_val = 7
min_E, min_A = find_minimums(K_val)

# Refine the search in a smaller window for better precision
def refine_search(K, E_center, A_center):
    E_range = np.linspace(E_center - 0.005, E_center + 0.005, 201)
    A_range = np.linspace(A_center - 0.01, A_center + 0.01, 201)
    
    min_E_found = float('inf')
    corresponding_A = -1

    for E in E_range:
        for A in A_range:
            if check_psd(E, A, K):
                if E < min_E_found:
                    min_E_found = E
                    corresponding_A = A
        # If we found any valid A for this E, it means this E is possible.
        # The tip of the allowed region is where the valid A range shrinks to a point.
        # We are looking for the lowest E that has at least one valid A.
        if min_E_found < float('inf'):
             # We can break here because we iterate E from low to high.
             break

    return min_E_found, corresponding_A

# Run the refined search if an initial point was found
if min_E != float('inf'):
    min_E, min_A = refine_search(K_val, min_E, min_A)


print(f"Minimal E: {min_E:.3f}")
print(f"Corresponding <x^2>: {min_A:.3f}")
