import math

def solve_diffie_hellman():
    """
    This function breaks a simple Diffie-Hellman key exchange by
    brute-forcing one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009
    g = 11

    # Publicly shared numbers from Alice and Bob
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    # We need to find either 'a' or 'b' to find the secret key.
    # Let's find Alice's secret number 'a'.
    # We are solving: 11^a mod 1009 = 297
    a = None
    for i in range(1, p):
        # Using Python's built-in pow(base, exp, mod) is efficient
        if pow(g, i, p) == alice_public:
            a = i
            break

    if a is not None:
        print(f"Found Alice's secret number 'a' by solving {g}^a mod {p} = {alice_public}.")
        print(f"Alice's secret number is a = {a}.")
        
        # Now, we can calculate the shared secret key using Bob's public number
        # and Alice's private number 'a'.
        # Secret Key = (bob_public)^a mod p
        secret_key = pow(bob_public, a, p)
        
        print("\nTo find the secret key, we calculate (Bob's public number)^a mod p.")
        print("The final equation is:")
        print(f"{bob_public}^{a} mod {p} = {secret_key}")
        
        print(f"\nThe secret key is {secret_key}.")
    else:
        print("Could not find the secret number 'a' through brute force.")

solve_diffie_hellman()
<<<330>>>