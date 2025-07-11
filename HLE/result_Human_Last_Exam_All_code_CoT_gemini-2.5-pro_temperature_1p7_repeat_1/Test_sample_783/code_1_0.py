def solve_diffie_hellman():
    """
    Breaks the Diffie-Hellman key exchange by brute-forcing the discrete logarithm.
    """
    p = 1009
    g = 11
    A = 297  # Alice's public number (g^a mod p)
    B = 944  # Bob's public number (g^b mod p)

    a = None
    # Step 1: Find Alice's secret number 'a' by brute force.
    # We are looking for 'a' such that g^a mod p = A
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break
    
    # Step 2: If 'a' is found, calculate the shared secret key.
    # The secret key S = B^a mod p
    if a is not None:
        secret_key = pow(B, a, p)
        
        # Print the final calculation as an equation
        print("The shared secret key can be computed using Bob's public number and Alice's private number.")
        print("Secret Key = (Bob's Public Number) ^ (Alice's Private Number) mod (Prime)")
        print(f"The equation is: {B} ^ {a} mod {p} = {secret_key}")
        print(f"\nThe secret key is: {secret_key}")
    else:
        print("Could not find Alice's secret number.")

solve_diffie_hellman()
<<<394>>>