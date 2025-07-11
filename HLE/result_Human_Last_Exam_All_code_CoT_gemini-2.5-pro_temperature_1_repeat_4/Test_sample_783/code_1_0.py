import sys

# Disable the creation of .pyc files
sys.dont_write_bytecode = True

def solve_diffie_hellman():
    """
    Solves for the secret key in a Diffie-Hellman exchange by brute-forcing
    the discrete logarithm problem for a small prime.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Generator

    # Publicly shared numbers from Alice and Bob
    A = 297  # Alice's public number (g^a mod p)
    B = 944  # Bob's public number (g^b mod p)

    print(f"Public information:")
    print(f"Prime p = {p}")
    print(f"Generator g = {g}")
    print(f"Alice's public number A = {A}")
    print(f"Bob's public number B = {B}\n")

    a = None
    b = None

    # Find Alice's secret number 'a' by solving g^a mod p = A
    print("Step 1: Finding Alice's secret number 'a'...")
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            print(f"Found Alice's secret number: a = {a}")
            break

    # Find Bob's secret number 'b' by solving g^b mod p = B
    print("Step 2: Finding Bob's secret number 'b'...")
    for i in range(1, p):
        if pow(g, i, p) == B:
            b = i
            print(f"Found Bob's secret number: b = {b}\n")
            break

    if a is None or b is None:
        print("Could not find secret numbers 'a' or 'b'.")
        return

    # Calculate the shared secret key
    print("Step 3: Calculating the shared secret key.")
    # Alice computes the secret key using Bob's public number and her secret number 'a'
    # S = B^a mod p
    secret_key = pow(B, a, p)

    print("Alice computes: S = (Bob's Public Key)^a mod p")
    print(f"The final equation is: {B}^{a} mod {p} = {secret_key}")

    # For verification, Bob can compute the same key using Alice's public number
    # and his secret number 'b'. The result must be the same.
    secret_key_from_b = pow(A, b, p)
    print("\nFor verification, Bob computes: S = (Alice's Public Key)^b mod p")
    print(f"The final equation is: {A}^{b} mod {p} = {secret_key_from_b}")

    if secret_key == secret_key_from_b:
        print(f"\nBoth calculations yield the same secret key: {secret_key}")
    else:
        print("\nError: The calculated keys do not match.")

if __name__ == '__main__':
    solve_diffie_hellman()
<<<324>>>