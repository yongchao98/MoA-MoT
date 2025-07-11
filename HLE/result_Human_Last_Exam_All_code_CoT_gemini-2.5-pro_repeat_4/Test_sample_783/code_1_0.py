def find_secret_key():
    """
    This function finds the secret key in a Diffie-Hellman exchange
    by first brute-forcing one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Generator

    # Publicly shared numbers
    alice_public = 297  # g^a mod p
    bob_public = 944    # g^b mod p

    # To find the shared secret, we first need to find either
    # Alice's secret 'a' or Bob's secret 'b'.
    # We can do this by brute-force since the prime 'p' is small.
    # Let's find Alice's secret 'a' by solving: 11^a mod 1009 = 297
    
    alice_secret = None
    for a in range(1, p):
        if pow(g, a, p) == alice_public:
            alice_secret = a
            break

    if alice_secret is None:
        print("Could not determine Alice's secret number.")
        return

    # Now that we have Alice's secret number 'a', we can compute the
    # shared secret key 'S' using Bob's public number: S = (g^b)^a mod p
    shared_secret = pow(bob_public, alice_secret, p)
    
    print("To find the shared secret key, we first find Alice's secret number 'a'.")
    print(f"By solving {g}^a mod {p} = {alice_public}, we find that a = {alice_secret}.")
    print("\nThen, we use 'a' and Bob's public number to find the shared secret key:")
    print(f"Equation: {bob_public} ^ {alice_secret} mod {p} = {shared_secret}")
    print(f"\nThe secret key is {shared_secret}")

find_secret_key()
<<<830>>>