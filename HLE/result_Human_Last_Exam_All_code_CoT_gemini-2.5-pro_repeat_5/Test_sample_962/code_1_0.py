import random
import hashlib

def main():
    """
    Main function to demonstrate the two-round threshold signature scheme.
    """
    # --- 1. SETUP ---
    # These are the parameters of our scheme.
    # n: The total number of participants (key share holders).
    # t: The threshold, i.e., the minimum number of participants required to sign.
    n = 5
    t = 3

    # We will perform all calculations in a finite field modulo a prime 'q'.
    # This prime should be large in a real-world scenario.
    q = 100003  # A sample prime number

    # The message to be signed.
    message = "This is a test message for the threshold signature scheme."

    print("--- 1. SETUP ---")
    print(f"Total participants (n): {n}")
    print(f"Threshold (t): {t}")
    print(f"Finite field prime (q): {q}")
    print(f"Message to sign: '{message}'")
    print("-" * 20 + "\n")


    # --- 2. KEY GENERATION (Performed by a Trusted Dealer) ---
    # The dealer creates a secret polynomial of degree t-1.
    # The master secret key 'sk' is the constant term, f(0).
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # where a_0 is the master secret key.
    poly = [random.randint(1, q - 1) for _ in range(t)]
    master_sk = poly[0]

    # The group public key 'pk' is derived from the master secret key.
    # In a real elliptic curve group, pk = sk * G (where G is the generator).
    # In our simplified model, the public key is just the number itself.
    master_pk = master_sk

    def evaluate_poly(x_coord):
        """Evaluates the polynomial f(x) at a given x coordinate."""
        y_coord = 0
        for i in range(t - 1, -1, -1):
            y_coord = (y_coord * x_coord + poly[i]) % q
        return y_coord

    # Generate 'n' secret shares by evaluating the polynomial at points x = 1, 2, ..., n
    secret_shares = {i: evaluate_poly(i) for i in range(1, n + 1)}

    print("--- 2. KEY GENERATION ---")
    print(f"Master Secret Key (sk = f(0)): {master_sk}")
    print(f"Group Public Key (pk): {master_pk}")
    print("\nSecret Shares distributed to participants:")
    for participant_id, share in secret_shares.items():
        print(f"  Participant {participant_id}: (x={participant_id}, y=f({participant_id})) = {share}")
    print("-" * 20 + "\n")


    # --- 3. SIGNING PROTOCOL ---
    # A subset of 't' participants will collaborate to sign the message.
    # Let's choose the first 't' participants for this demonstration.
    signer_ids = list(range(1, t + 1))
    
    print("--- 3. SIGNING PROTOCOL ---")
    print(f"Signing Participants (S): {signer_ids}")

    # --- ROUND 1: Commitment ---
    # Each participant 'i' in the signing group generates a secret random nonce 'k_i'
    # and computes a public commitment 'R_i'. In our simplified model, R_i = k_i.
    nonces = {}
    commitments = {}
    print("\n--- Round 1: Commitment ---")
    for i in signer_ids:
        k_i = random.randint(1, q - 1)
        R_i = k_i
        nonces[i] = k_i
        commitments[i] = R_i
        print(f"  Participant {i}:")
        print(f"    Secret Nonce (k_{i}): {k_i}")
        print(f"    Public Commitment (R_{i}): {R_i}")

    # --- ROUND 2: Signature Share Generation ---
    # All participants receive the commitments and proceed.
    print("\n--- Round 2: Signature Share Generation ---")
    
    # a) All participants compute the group commitment 'R'.
    group_commitment_R = sum(commitments.values()) % q
    print(f"Group Commitment (R = sum(R_i)): {group_commitment_R}")

    # b) All participants compute the challenge 'c'.
    # c = H(pk, R, m) where H is a cryptographic hash function.
    hasher = hashlib.sha256()
    hasher.update(str(master_pk).encode())
    hasher.update(str(group_commitment_R).encode())
    hasher.update(message.encode())
    challenge_c = int(hasher.hexdigest(), 16) % q
    print(f"Challenge (c = H(pk, R, m)): {challenge_c}")

    def mod_inverse(n, prime):
        """Computes the modular inverse of n modulo prime using Fermat's Little Theorem."""
        return pow(n, prime - 2, prime)

    def get_lagrange_coeff(i, S):
        """
        Calculates the Lagrange basis polynomial l_i(0) for participant i in set S.
        l_i(0) = product(j / (j - i) for j in S if j != i)
        """
        num = 1
        den = 1
        for j in S:
            if i == j:
                continue
            num = (num * j) % q
            den = (den * (j - i)) % q
        return (num * mod_inverse(den, q)) % q

    # c) Each participant computes their signature share 'z_i'.
    # z_i = k_i + c * lambda_i * s_i
    # where s_i is the participant's secret share and lambda_i is their Lagrange coefficient.
    sig_shares = {}
    print("\nSignature Shares:")
    for i in signer_ids:
        s_i = secret_shares[i]
        k_i = nonces[i]
        lambda_i = get_lagrange_coeff(i, signer_ids)
        z_i = (k_i + challenge_c * lambda_i * s_i) % q
        sig_shares[i] = z_i
        print(f"  Participant {i}:")
        print(f"    Lagrange Coeff (lambda_{i}): {lambda_i}")
        print(f"    Signature Share (z_{i} = k_{i} + c*lambda_{i}*s_{i}): {z_i}")
    print("-" * 20 + "\n")


    # --- 4. AGGREGATION AND VERIFICATION ---
    print("--- 4. AGGREGATION & VERIFICATION ---")
    
    # a) The Aggregator sums the shares to create the final signature component 'z'.
    final_z = sum(sig_shares.values()) % q
    
    # The final signature is the pair (R, z).
    final_signature = (group_commitment_R, final_z)
    print(f"Final Aggregated Signature (R, z): {final_signature}")

    # b) Verification
    # To verify, anyone can check if z == R + c * pk (mod q).
    # First, the verifier must compute the challenge 'c' exactly as the signers did.
    verifier_hasher = hashlib.sha256()
    verifier_hasher.update(str(master_pk).encode())
    verifier_hasher.update(str(final_signature[0]).encode()) # Use R from the signature
    verifier_hasher.update(message.encode())
    verifier_challenge_c = int(verifier_hasher.hexdigest(), 16) % q

    # The left side of the verification equation is 'z'.
    lhs = final_signature[1]
    
    # The right side is 'R + c * pk'.
    rhs = (final_signature[0] + verifier_challenge_c * master_pk) % q

    print("\nVerification Equation: z == R + c * pk (mod q)")
    print(f"  z (Left Side): {lhs}")
    print(f"  R + c*pk (Right Side): {rhs}")

    is_valid = (lhs == rhs)
    print(f"\nSignature validity: {is_valid}")
    if is_valid:
        print("✅ The signature is valid!")
    else:
        print("❌ The signature is invalid!")

if __name__ == '__main__':
    main()
