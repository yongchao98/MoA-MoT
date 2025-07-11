import numpy as np

def get_sun_generators(N):
    """
    Generates a basis of N^2 - 1 traceless, Hermitian generators for SU(N),
    normalized such that Tr(T_a T_b) = (1/2) * delta_ab.
    This uses a generalization of the Gell-Mann matrices.
    """
    if not isinstance(N, int) or N < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    num_generators = N**2 - 1
    generators = []

    # Type 1: Symmetric matrices
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = 1
            mat[k, j] = 1
            # Normalize to Tr(T^2) = 1/2. Tr(mat^2) = 2.
            # So T = mat / sqrt(2*Tr(mat^2)) = mat / 2
            generators.append(mat / 2.0)

    # Type 2: Anti-symmetric matrices
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = -1j
            mat[k, j] = 1j
            # Normalize to Tr(T^2) = 1/2. Tr(mat^2) = 2.
            # So T = mat / 2
            generators.append(mat / 2.0)

    # Type 3: Diagonal matrices
    for l in range(1, N):
        mat = np.zeros((N, N), dtype=complex)
        diag_elements = [1.0] * l + [-l] + [0.0] * (N - l - 1)
        # Normalization factor for the diagonal matrix part to have Tr(D^2)=2
        norm = np.sqrt(2 * l * (l + 1))
        # Overall normalization factor to have Tr(T^2)=1/2
        # Normalization factor is 1/2 * np.sqrt(2/(l*(l+1))) is wrong.
        # Lambda must have Tr(L^2) = 2. T = L/2.
        # D_l = 1/sqrt(l(l+1)/2) * diag(...) so Tr(D_l^2)=2
        # Our diag elements give trace of l*(1^2)+(-l)^2 = l+l^2=l(l+1)
        # So we must divide by sqrt(l(l+1)) to get Tr(D^2)=2.
        # Hence T = D/2
        final_norm = 1.0 / np.sqrt(l * (l + 1.0)) * (1.0/np.sqrt(2))
        # Let's verify Tr(T^2)
        # T_diag = np.diag(diag_elements) * (1/sqrt(2*l*(l+1)))
        # T_diag^2 has diag elems [1, 1, ... l^2] / (2*l*(l+1))
        # Trace is (l + l^2) / (2*l*(l+1)) = 1/2. Correct.
        factor = 1.0 / np.sqrt(2 * l * (l + 1))
        for m in range(l):
            mat[m,m] = factor
        mat[l,l] = -l * factor

        generators.append(mat)
    
    # Assert correct number of generators were created
    assert len(generators) == num_generators
    
    return generators

def count_distinct_d_values(N, precision=8):
    """
    Calculates the number of distinct non-zero symmetric structure constants d_ijk for SU(N).
    
    Args:
    N (int): The dimension of the special unitary group.
    precision (int): The number of decimal places to round to for uniqueness check.
    
    Returns:
    int: The number of distinct non-zero d_ijk values.
    """
    # Get the SU(N) generators
    T = get_sun_generators(N)
    num_gen = len(T)
    
    # Use a set to store unique d_ijk values
    d_values = set()
    
    # Set a tolerance for checking if a value is non-zero
    tolerance = 10**(-precision-1)
    
    # Iterate over all combinations of i, j, k with i <= j <= k
    for i in range(num_gen):
        for j in range(i, num_gen):
            for k in range(j, num_gen):
                # Formula for d_ijk derived from the anti-commutation relation
                # d_ijk = Tr(T_i T_j T_k) + Tr(T_j T_i T_k)
                # Since T_a are Hermitian, the result must be real. We take the real part
                # to discard any small imaginary noise from floating point inaccuracies.
                d_val = np.trace(T[i] @ T[j] @ T[k]) + np.trace(T[j] @ T[i] @ T[k])
                d_val = d_val.real
                
                # Check if the value is significantly different from zero
                if abs(d_val) > tolerance:
                    # Round to handle floating point variations of the same number
                    rounded_val = round(d_val, precision)
                    d_values.add(rounded_val)

    return len(d_values)

if __name__ == '__main__':
    # You can change N to the desired value.
    # N=2 -> 0 values
    # N=3 -> 4 values
    # N=4 -> 9 values
    N = 4 
    
    print(f"For SU({N}), the generators T_i are defined to satisfy:")
    print(f"{{T_i, T_j}} = (1/{N}) * delta_ij * I + 2 * sum_k(d_ijk * T_k)")
    print(f"and Tr(T_i * T_j) = 0.5 * delta_ij.\n")

    num_values = count_distinct_d_values(N)

    print(f"The number of different numerical values the non-zero d_ijk take for SU({N}) is: {num_values}")