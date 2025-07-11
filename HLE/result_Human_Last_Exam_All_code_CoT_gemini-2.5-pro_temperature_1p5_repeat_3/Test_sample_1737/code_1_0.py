import numpy as np

def solve_su_n_structure_constants():
    """
    Calculates the number of unique non-zero values for the totally symmetric
    structure constants d_ijk of SU(N).
    """
    try:
        N = int(input("Enter the value of N for SU(N): "))
        if N < 2:
            print("N must be an integer greater than or equal to 2.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The dimension of the Lie algebra su(n) is n^2 - 1.
    dim = N * N - 1
    if dim == 0:
        print("For N=1, the algebra is trivial and there are no structure constants.")
        print("Number of different numerical values: 0")
        return

    # Generate the generalized Gell-Mann matrices (lambda_a)
    # These matrices are hermitian, traceless, and satisfy Tr(L_a L_b) = 2 * delta_ab
    lambdas = []

    # 1. Symmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            u_jk = np.zeros((N, N), dtype=complex)
            u_jk[j, k] = 1.0
            u_jk[k, j] = 1.0
            lambdas.append(u_jk)

    # 2. Anti-symmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            v_jk = np.zeros((N, N), dtype=complex)
            v_jk[j, k] = -1.0j
            v_jk[k, j] = 1.0j
            lambdas.append(v_jk)

    # 3. Diagonal matrices
    for l in range(1, N):
        w_l = np.zeros((N, N), dtype=complex)
        for m in range(l):
            w_l[m, m] = 1.0
        w_l[l, l] = -float(l)
        w_l *= np.sqrt(2.0 / (l * (l + 1)))
        lambdas.append(w_l)

    # The SU(N) generators T_a are normalized as Tr(T_a T_b) = 0.5 * delta_ab.
    # This is achieved by T_a = lambda_a / 2.
    generators = [l / 2.0 for l in lambdas]

    # Calculate d_ijk = 4 * Re(Tr(T_i T_j T_k))
    # We use a set to store unique values.
    d_values = set()
    tolerance = 1e-9 # Tolerance for floating point zero
    precision = 8    # Rounding precision

    # Iterate i <= j <= k due to total symmetry of d_ijk
    for i in range(dim):
        for j in range(i, dim):
            for k in range(j, dim):
                # Matrix multiplication and trace
                trace_val = np.trace(generators[i] @ generators[j] @ generators[k])
                
                # The d_ijk constant
                d_val = 4 * np.real(trace_val)

                # Round to handle floating point inaccuracies
                d_val_rounded = round(d_val, precision)
                
                # Add to set if non-zero
                if abs(d_val_rounded) > tolerance:
                    d_values.add(d_val_rounded)
    
    sorted_values = sorted(list(d_values))
    print(f"\nFor SU({N}), the unique non-zero d_ijk values are:")
    # "output each number in the final equation" is interpreted as printing the found values
    for val in sorted_values:
        print(val)
        
    print(f"\nThe number of different numerical values for the non-zero d_ijk is:")
    print(len(d_values))

if __name__ == '__main__':
    solve_su_n_structure_constants()