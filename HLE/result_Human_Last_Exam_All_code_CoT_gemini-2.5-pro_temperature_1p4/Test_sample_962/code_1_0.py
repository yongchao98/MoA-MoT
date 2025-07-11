import random
import hashlib

def design_threshold_signature_scheme(t, n, message):
    """
    Simulates a 2-round t-out-of-n threshold signature scheme.

    Args:
        t (int): The threshold of signers required.
        n (int): The total number of parties in the scheme.
        message (str): The message to be signed.
    """
    if t > n:
        print("Threshold t cannot be greater than the number of parties n.")
        return

    print(f"Designing a {t}-out-of-{n} threshold signature scheme.\n")

    # --- 1. SETUP PHASE ---
    # In a real system, this is a complex Distributed Key Generation (DKG) protocol.
    # We simulate it locally for demonstration.

    # Choose a large prime for our finite field (in ECC, this would be the curve order)
    q = 1009 # A small prime for demonstration

    # Generate a secret polynomial of degree t-1
    # f(x) = a_0 + a_1*x + a_2*x^2 + ... + a_{t-1}*x^{t-1}
    # The master secret key is s = f(0) = a_0
    coefficients = [random.randint(1, q - 1) for _ in range(t - 1)]
    s = random.randint(1, q - 1) # This is a_0, the master secret key
    coefficients.insert(0, s)

    def f(x):
        """The secret polynomial."""
        y = 0
        for i, coeff in enumerate(coefficients):
            y += coeff * (x ** i)
        return y % q

    # Generate secret shares for n parties
    # Party i is given the secret share s_i = f(i)
    secret_shares = {i: f(i) for i in range(1, n + 1)}
    print("--- 1. Key Generation (Setup) ---")
    print(f"Chosen prime q: {q}")
    print(f"Master secret key (s = f(0)): {s}")
    # print(f"Secret shares distributed (s_i = f(i)): {secret_shares}")
    print(f"{n} secret shares have been generated and distributed securely.")
    print("-" * 20)

    # --- 2. SIGNING PHASE ---
    # A subset of t parties will now sign the message.
    signer_indices = list(range(1, t + 1))
    print(f"\n--- 2. Signing Protocol ---")
    print(f"Message to sign: '{message}'")
    print(f"Signers (first {t} parties): {signer_indices}\n")

    # -- ROUND 1: Commitment --
    print("-- Round 1: Commitment --")
    nonces = {}
    commitments = {}
    for i in signer_indices:
        # Each signer generates a secret random nonce r_i
        r_i = random.randint(1, q - 1)
        nonces[i] = r_i
        # and broadcasts a commitment R_i. In real ECC, R_i = r_i * G.
        # Here, we just use the nonce itself for simplicity.
        commitments[i] = r_i
        print(f"Party {i} generates a secret nonce and broadcasts commitment R_{i} = {commitments[i]}")

    # -- ROUND 2: Response --
    print("\n-- Round 2: Response --")
    # All parties received the commitments. They compute an aggregate commitment R.
    R = sum(commitments.values()) % q
    print(f"All signers compute the aggregate commitment R = sum(R_i) = {R}")

    # All parties compute a challenge 'c' based on R and the message.
    # A cryptographic hash function must be used.
    hash_input = f"{message}{R}".encode('utf-8')
    c = int(hashlib.sha256(hash_input).hexdigest(), 16) % q
    print(f"All signers compute the challenge c = H(message, R) = {c}")

    # Helper function for modular inverse, needed for Lagrange coefficients.
    def modInverse(k, mod):
        return pow(k, -1, mod)

    # Helper for Lagrange interpolation at x=0.
    def get_lagrange_coeff(i, S):
        num = 1
        den = 1
        for j in S:
            if i != j:
                num = (num * j)
                den = (den * (j - i))
        return (num * modInverse(den, q)) % q

    # Each signer computes their partial signature z_i
    partial_signatures = {}
    print("\nEach signer computes and broadcasts their partial signature z_i:")
    for i in signer_indices:
        s_i = secret_shares[i]
        r_i = nonces[i]
        # Calculate this signer's Lagrange coefficient for x=0
        lambda_i = get_lagrange_coeff(i, signer_indices)
        # Calculate the partial signature
        z_i = (r_i + c * lambda_i * s_i) % q
        partial_signatures[i] = z_i
        print(f"Party {i}: z_{i} = (r_{i} + c * lambda_{i} * s_{i}) mod q = {z_i}")

    # --- 3. AGGREGATION & VERIFICATION ---
    print("\n--- 3. Aggregation and Verification ---")
    # An aggregator combines the partial signatures to get the final signature (R, Z)
    Z = sum(partial_signatures.values()) % q
    signature = (R, Z)
    print(f"Aggregated Signature (R, Z): {signature}")

    # A verifier checks the signature using the master public key `s`.
    # In a real system, they would use the public key P = s*G.
    print("\n--- Final Equation Check: Z = (R + c*s) mod q ---")
    
    # Calculate the right-hand side of the verification equation
    RHS = (R + c * s) % q

    print(f"Master Secret (s): {s}")
    print(f"Aggregate Nonce (R): {R}")
    print(f"Challenge (c): {c}")
    print(f"Aggregate Response (Z): {Z}")
    
    print("\nEquation Check:")
    print(f"LHS (Z) = {Z}")
    print(f"RHS (R + c*s) = {RHS}")

    is_valid = (Z == RHS)
    print(f"\nSignature is valid: {is_valid}")
    if not is_valid:
        print("Verification failed! Something is wrong.")


if __name__ == '__main__':
    # --- Parameters ---
    N_PARTIES = 5 # Total number of parties
    THRESHOLD = 3 # Number of parties needed to sign
    MESSAGE_TO_SIGN = "Alice pays Bob 10 BTC"

    design_threshold_signature_scheme(THRESHOLD, N_PARTIES, MESSAGE_TO_SIGN)