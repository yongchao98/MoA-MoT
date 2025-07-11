def solve_diffie_hellman():
    """
    Solves for the shared secret in a Diffie-Hellman key exchange
    by brute-forcing one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009
    g = 11

    # Public keys exchanged
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    # We will find Alice's secret key 'a' by brute force.
    # We can iterate from 1 up to p-1.
    alice_secret = None
    for a in range(1, p):
        # Calculate g^a mod p
        if pow(g, a, p) == alice_public:
            alice_secret = a
            break

    if alice_secret is not None:
        # Now that we have Alice's secret key 'a', we can calculate the shared secret.
        # The shared secret S = (Bob's public key)^a mod p
        shared_secret = pow(bob_public, alice_secret, p)

        print("The prime number (p) is: 1009")
        print("The base number (g) is: 11")
        print("Alice's public key (A) is: 297")
        print("Bob's public key (B) is: 944")
        print(f"Found Alice's secret key (a): {alice_secret}")
        print("\nCalculating the shared secret key...")
        print("Formula: S = (Bob's Public Key) ^ (Alice's Secret Key) mod p")
        # Output the final equation with all the numbers
        print(f"Final Equation: {bob_public}^{alice_secret} mod {p} = {shared_secret}")
        print(f"\nThe secret key is: {shared_secret}")
    else:
        print("Could not find Alice's secret key through brute force.")

solve_diffie_hellman()
<<<16>>>