def compute_generalized_markov_number_cf():
    """
    Computes the generalized Markov number m_{4/7} and its continued fraction.
    """
    p, q = 4, 7

    # Step 1: Compute the continued fraction coefficients for p/q = [0; a1, a2, ...].
    # We compute the coefficients for q/p = [a1; a2, ...].
    coeffs = []
    num, den = q, p
    while den > 0:
        a = num // den
        coeffs.append(a)
        num, den = den, num % den
    
    # Step 2: Construct the Markov matrix M.
    # M = A_n^a_n * ... * A_1^a_1
    # A_k = R if k is odd, L if k is even.
    # For R^k, the matrix is [[1, k], [0, 1]].
    # For L^k, the matrix is [[1, 0], [k, 1]].

    def mat_mul(A, B):
        """Multiplies two 2x2 matrices."""
        C = [[0, 0], [0, 0]]
        C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
        C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
        C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
        C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
        return C

    n = len(coeffs)
    
    # Initialize M with the rightmost term in the product: A_1^a_1
    # k=1 is odd, so A_1 = R.
    a1 = coeffs[0]
    M = [[1, a1], [0, 1]]  # R^a1

    # Multiply by the other terms from right to left (i.e., k=2 to n)
    for i in range(1, n):
        k = i + 1
        a_k = coeffs[i]
        
        if k % 2 == 0:  # k is even, A_k = L
            A_k_power = [[1, 0], [a_k, 1]]  # L^a_k
        else:  # k is odd, A_k = R
            A_k_power = [[1, a_k], [0, 1]]  # R^a_k
            
        M = mat_mul(A_k_power, M)

    # Step 3: Calculate the generalized Markov number.
    trace_M = M[0][0] + M[1][1]
    markov_number = trace_M / 3
    
    # The result should be an integer.
    markov_number_int = int(round(markov_number))

    # Step 4: Compute the continued fraction of the resulting number.
    # For an integer N, the continued fraction is simply [N].
    final_cf = [markov_number_int]

    # Step 5: Print the final equation.
    cf_str = '[' + '; '.join(map(str, final_cf)) + ']'
    print(f"The generalized Markov number m_{p}/{q} is {markov_number_int}.")
    print("Its continued fraction is given by the equation:")
    print(f"{markov_number_int} = {cf_str}")

compute_generalized_markov_number_cf()