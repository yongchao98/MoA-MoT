import numpy as np

def solve_bootstrap():
    """
    Solves for the ground state energy E and <x^2> of the V(x) = x^2 + x^4 potential
    using the quantum bootstrap method with K=7.
    """
    # K=7 implies we need moments up to <x^(7+7)> = <x^14>
    K = 7
    MAX_MOMENT_IDX = 2 * K

    def calculate_moments(E, v2):
        """
        Calculates the moments <x^k> for a given E and v2 = <x^2>.
        The potential is V(x) = x^2 + x^4.
        """
        if v2 < 0:
            return None
        
        # Initialize moments dictionary
        # v_n = <x^n>
        v = {0: 1.0, 2: v2}
        for i in range(1, MAX_MOMENT_IDX + 2, 2):
            v[i] = 0.0

        # The core recursion relation is derived from <[H, x^t p]> = 0
        # v_{t+3} = (1/(4t+8)) * [4tE*v_{t-1} - (4t+4)*v_{t+1} + t(t-1)(t-2)*v_{t-3}]
        # We use this for odd t, starting from t=1, to find all even moments.

        # Calculate v[4] using t=1
        v[4] = (4*1*E*v[0] - (4*1+4)*v[2] + 0) / (4*1 + 8)
        v[4] = (4*E - 8*v2) / 12.0

        # Calculate remaining even moments up to v[14] using the recursion for t = 3, 5, ...
        for t in range(3, MAX_MOMENT_IDX - 1, 2):
            idx_p3 = t + 3
            if idx_p3 > MAX_MOMENT_IDX:
                break
            
            v_t_minus_3 = v.get(t - 3, 0)
            v_t_minus_1 = v.get(t - 1, 0)
            v_t_plus_1 = v.get(t + 1, 0)

            numerator = (4*t*E*v_t_minus_1 - (4*t+4)*v_t_plus_1 + t*(t-1)*(t-2)*v_t_minus_3)
            denominator = 4*t + 8
            
            if denominator == 0: return None # Avoid division by zero
            v[idx_p3] = numerator / denominator
        
        return v

    def check_positivity(moments):
        """
        Checks if the moment matrices M_even and M_odd are positive semidefinite.
        """
        if moments is None:
            return False
            
        # For K=7, the basis for even operators is {1, x^2, x^4, x^6}
        # and for odd operators is {x, x^3, x^5, x^7}.
        # This results in two 4x4 matrices.
        dim = (K // 2) + 1

        # Build M_even, M_ij = <x^(2i) x^(2j)> = v_{2i+2j}
        M_even = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(dim):
                M_even[i, j] = moments[2*i + 2*j]
                
        # Build M_odd, M_ij = <x^(2i+1) x^(2j+1)> = v_{2i+2j+2}
        M_odd = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(dim):
                M_odd[i, j] = moments[2*i + 2*j + 2]
        
        try:
            # Check if all eigenvalues are non-negative (within a small tolerance)
            eigvals_even = np.linalg.eigvalsh(M_even)
            if np.any(eigvals_even < -1e-9):
                return False

            eigvals_odd = np.linalg.eigvalsh(M_odd)
            if np.any(eigvals_odd < -1e-9):
                return False
        except np.linalg.LinAlgError:
            # Catches errors from non-finite matrix elements (NaN, Inf)
            return False
            
        return True

    def find_allowed_v2(E):
        """
        For a given E, scan for a v2 that satisfies the positivity constraints.
        Returns a valid v2 if found, otherwise None.
        """
        # A necessary (but not sufficient) condition comes from the 2x2 M_even matrix,
        # which sets an upper bound on v2 for a given E.
        # 3*v2^2 + 2*v2 - E <= 0
        if 4 + 12 * E < 0:
            return None
        v2_max = (-2 + np.sqrt(4 + 12 * E)) / 6.0
        
        # Scan v2 in its allowed range to find if a valid solution exists.
        for v2 in np.linspace(0.001, v2_max, 250):
            moments = calculate_moments(E, v2)
            if check_positivity(moments):
                return v2 
        return None

    # --- Main Search Logic ---
    # Binary search for the minimal allowed energy E
    E_low, E_high = 1.0, 2.0
    min_E = -1
    
    for _ in range(25): # 25 iterations provide high precision
        E_mid = (E_low + E_high) / 2.0
        if find_allowed_v2(E_mid) is not None:
            # E_mid is a possible energy, so we try to find an even lower one
            min_E = E_mid
            E_high = E_mid
        else:
            # E_mid is not a possible energy, so the minimum must be higher
            E_low = E_mid
    
    # Refined search for the corresponding <x^2> at the minimal energy
    final_E = min_E
    v2_max = (-2 + np.sqrt(4 + 12 * final_E)) / 6.0
    final_v2 = -1
    
    # At the minimal E, the allowed range for v2 should be very small (a single point).
    # We do a fine-grained search to locate it.
    for v2 in np.linspace(0.001, v2_max, 5001):
        moments = calculate_moments(final_E, v2)
        if check_positivity(moments):
            final_v2 = v2
            break

    # Print the final result
    print(f"Minimal E: {final_E:.3f}")
    print(f"Minimal <x^2>: {final_v2:.3f}")

    return final_E, final_v2

# --- Execute the Solver ---
E_result, v2_result = solve_bootstrap()
answer_string = f"E = {E_result:.3f}, <x^2> = {v2_result:.3f}"
print(f"\n<<<{answer_string}>>>")