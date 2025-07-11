import numpy as np

def solve_bootstrap_qm():
    """
    Solves for the ground state energy (E) and <x^2> of the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """
    K = 7
    # Use memoization (caching) to store results of expensive calculations
    # and avoid re-computing them with the same inputs.
    memo_moments = {}
    memo_psd_check = {}

    def calculate_moments(E, m2):
        """
        Calculates even moments <x^{2k}> for k from 0 to K using the recursion relation.
        <x^0> and <x^2> are given as inputs.
        """
        # Return cached result if available
        if (E, m2) in memo_moments:
            return memo_moments[(E, m2)]

        # Array to store even moments: m_even[k] = <x^{2k}>
        # For K=7, we need moments up to <x^{2*K}> = <x^{14}>. The array size is K+1.
        m_even = np.zeros(K + 1)
        m_even[0] = 1.0  # <x^0> = 1
        m_even[1] = m2   # <x^2> is a test parameter

        # The recursion from Step 3 for V(x) = x^2 + x^4 simplifies to:
        # <x^{2k+2}> = 1/(4(2k+1)) * [ 4(2k-1)E*<x^{2k-2}> - 8k*<x^{2k}> + (2k-1)(2k-2)(2k-3)*<x^{2k-4}> ]
        # This formula is valid for k >= 1.
        for k in range(1, K): # Loop from k=1 to 6 to calculate m_even[2] through m_even[7]
            m2k = m_even[k]
            m2k_minus_2 = m_even[k-1]
            m2k_minus_4 = m_even[k - 2] if k >= 2 else 0.0

            term3 = (2 * k - 1) * (2 * k - 2) * (2 * k - 3) * m2k_minus_4
            numerator = 4 * (2 * k - 1) * E * m2k_minus_2 - 8 * k * m2k + term3
            denominator = 4 * (2 * k + 1)
            
            m_even[k + 1] = numerator / denominator
        
        # Cache and return the calculated moments
        memo_moments[(E, m2)] = m_even
        return m_even

    def check_psd(E, m2):
        """
        Checks if the moment matrices M_even and M_odd are positive semi-definite (PSD).
        A matrix is PSD if all its eigenvalues are non-negative.
        """
        # Return cached result if available
        if (E, m2) in memo_psd_check:
            return memo_psd_check[(E, m2)]

        m_even_values = calculate_moments(E, m2)

        # If moments calculation led to invalid numbers, it's not a valid point.
        if np.any(np.isinf(m_even_values)) or np.any(np.isnan(m_even_values)):
            return False
            
        # For K=7, the basis operators are {1, x, x^2, ..., x^7}.
        # Due to symmetry, this splits into even and odd operator bases.
        # Even basis: {1, x^2, x^4, x^6}, giving a 4x4 matrix M_even.
        # Odd basis:  {x, x^3, x^5, x^7}, giving a 4x4 matrix M_odd.
        
        M_even = np.zeros((4, 4))
        M_odd = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                # M_even_(i,j) = <x^{2i} * x^{2j}> = <x^{2(i+j)}>
                M_even[i, j] = m_even_values[i + j]
                # M_odd_(i,j) = <x^{2i+1} * x^{2j+1}> = <x^{2(i+j+1)}>
                M_odd[i, j] = m_even_values[i + j + 1]

        try:
            # Use np.linalg.eigvalsh as matrices are symmetric, it's faster.
            eig_even = np.linalg.eigvalsh(M_even)
            eig_odd = np.linalg.eigvalsh(M_odd)
        except np.linalg.LinAlgError:
            # Matrix was not well-behaved (e.g., singular)
            return False

        # A point is valid if all eigenvalues of both matrices are non-negative.
        # We use a small negative tolerance to account for floating-point inaccuracies.
        tolerance = -1e-9
        result = np.min(eig_even) >= tolerance and np.min(eig_odd) >= tolerance
        memo_psd_check[(E, m2)] = result
        return result

    def has_any_valid_x2(E):
        """Scans a range of <x^2> values to see if any are valid for a given E."""
        x2_scan_range = np.linspace(0.2, 0.6, 81) # Scan plausible <x^2> values
        for x2_val in x2_scan_range:
            if check_psd(E, x2_val):
                return True # Found a valid <x^2>, so this E is possible
        return False # No valid <x^2> found for this E

    # --- Main Search Logic ---

    # 1. Binary search for the minimum possible energy, E.
    E_low, E_high = 1.0, 1.5  # Plausible range for the ground state energy
    for _ in range(25): # 25 iterations for high precision
        E_mid = (E_low + E_high) / 2
        if has_any_valid_x2(E_mid):
            E_high = E_mid  # A solution exists, so try for an even lower E
        else:
            E_low = E_mid   # No solution, E is too low
    min_E = E_high

    # 2. For this minimal E, find the corresponding value of <x^2>.
    # At the tip of the allowed region, there should be a unique <x^2>.
    # We can find it by scanning for the first valid value.
    min_x2 = -1.0
    x2_fine_scan = np.linspace(0.2, 0.6, 4001) # Very fine scan for accuracy
    for x2_val in x2_fine_scan:
        if check_psd(min_E, x2_val):
            min_x2 = x2_val
            break # Found the minimal <x^2> for our minimal E
    
    print(f"Minimal Energy (E): {min_E:.3f}")
    print(f"Corresponding <x^2>: {min_x2:.3f}")
    
    return min_E, min_x2

# Execute the solver
minimal_E, minimal_x2 = solve_bootstrap_qm()
# The final values formatted for the specific request
print(f"<<<E = {minimal_E:.3f}, <x^2> = {minimal_x2:.3f}>>>")