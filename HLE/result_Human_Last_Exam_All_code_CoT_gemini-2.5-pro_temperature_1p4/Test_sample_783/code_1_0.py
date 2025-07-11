def solve_diffie_hellman():
    """
    Breaks the Diffie-Hellman key exchange for the given parameters by
    brute-forcing one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Publicly shared numbers
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    print("--- Given Information ---")
    print(f"Prime modulus (p): {p}")
    print(f"Base (g): {g}")
    print(f"Alice's public number (A): {A}")
    print(f"Bob's public number (B): {B}")
    print("-" * 25)

    # To break the code, we need to find Alice's secret number 'a' or Bob's 'b'.
    # We can find 'a' by solving g^a mod p = A.
    # Since p is small, we can iterate through possible values of 'a'.
    a = None
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break

    if a is not None:
        print(f"Step 1: Found Alice's secret number 'a' by solving {g}^a mod {p} == {A}.")
        print(f"Alice's secret number (a) is: {a}")

        # Now we can compute the shared secret key 'S' using B^a mod p
        secret_key = pow(B, a, p)

        print("\nStep 2: Calculate the shared secret key (S) using B^a mod p.")
        print("The final equation to get the secret key is:")
        print(f"{B} ^ {a} mod {p} = {secret_key}")

        print("\n--- Result ---")
        print(f"The secret key is: {secret_key}")
    else:
        print("Could not find the secret number 'a' through brute force.")

solve_diffie_hellman()
<<<419>>>