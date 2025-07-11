def solve_diffie_hellman():
    """
    Breaks the Diffie-Hellman key exchange by finding the secret key.
    
    This is done by solving the discrete logarithm problem for a small prime
    using a brute-force search.
    """
    # Publicly known numbers
    p = 1009  # Prime number
    g = 11    # Generator
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # We need to find Alice's secret number 'a'.
    # We can do this by iterating through possible values of 'a' until we
    # find one that satisfies the equation g^a mod p = A.
    a = None
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break
    
    if a is None:
        print("Could not find Alice's secret number 'a'.")
        return

    # Now that we have found Alice's secret number 'a', we can calculate
    # the shared secret key 'S' using Bob's public number 'B'.
    # The formula is S = B^a mod p.
    secret_key = pow(B, a, p)

    print("Step 1: Find Alice's secret number 'a' by solving g^a mod p = A.")
    print(f"Solving the equation: {g}^a mod {p} = {A}")
    print(f"Found Alice's secret number: a = {a}\n")
    
    print("Step 2: Compute the shared secret key 'S' using S = B^a mod p.")
    print("The final equation is:")
    print(f"S = {B}^{a} mod {p}")
    print(f"The shared secret key is: {secret_key}")

solve_diffie_hellman()
<<<694>>>