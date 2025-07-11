import random
import hashlib
from math import prod

def design_threshold_signature_scheme():
    """
    Simulates a 2-round t-out-of-n threshold signature scheme based on FROST.
    """

    # Helper function for modular inverse using Python's built-in power function
    # pow(b, -1, q) is equivalent to b^(q-2) mod q for a prime q
    def mod_inverse(b, q):
        return pow(b, -1, q)

    # Helper for a consistent hash function (SHA-256)
    def H(*args):
        """A simple hash function wrapper."""
        s = "".join(map(str, args))
        return int(hashlib.sha256(s.encode()).hexdigest(), 16)

    # --- Elliptic Curve Simulation (using modulo arithmetic) ---
    # We use a large prime number (the order of the secp256k1 curve) as our field
    Q = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    # G is a "generator point" on our simulated curve
    G = 2

    # --- 1. Protocol Parameters & Key Generation (Simulated DKG) ---
    N = 5  # Total number of parties in the group
    T = 3  # The threshold of signers required

    print(f"--- Step 1: Setup ({T}-out-of-{N} Scheme) ---")
    
    # In a real DKG, parties collaboratively generate the secret key without a trusted dealer.
    # For this simulation, we'll create a secret polynomial of degree T-1.
    # f(z) = a_0 + a_1*z + ... + a_{T-1}*z^{T-1}
    # The master secret key is x = f(0) = a_0
    coeffs = [random.randrange(1, Q) for _ in range(T)]
    secret_key_x = coeffs[0]

    def secret_polynomial(z):
        """Evaluates the secret polynomial at a given point z."""
        return sum(coeffs[i] * (z ** i) for i in range(T)) % Q

    # Generate private key shares for each of the N parties. Party `i` gets share x_i = f(i).
    # The party identifiers are 1, 2, ..., N.
    private_key_shares = {i: secret_polynomial(i) for i in range(1, N + 1)}
    
    # The single public key for the group is P = x * G
    public_key_P = (secret_key_x * G) % Q
    print(f"Simulated DKG complete. Group public key P has been generated.")
    print("-" * 20)

    # --- 2. Signing Protocol ---
    # A subset of T parties will sign a message.
    # We randomly select T participants from the N parties.
    all_party_ids = list(range(1, N + 1))
    random.shuffle(all_party_ids)
    signers_S = set(all_party_ids[:T])
    
    message_m = "Design a secure two-round threshold signature scheme."
    print(f"--- Step 2: Signing Protocol Begins ---")
    print(f"Message to sign: '{message_m}'")
    print(f"Participating signers (party IDs): {sorted(list(signers_S))}")

    # --- ROUND 1: COMMITMENT ---
    print("\n--- Signing Round 1: Commitment ---")
    signer_nonces_k = {}
    signer_public_nonces_R = {}
    commitments_C = {}

    for i in signers_S:
        # Each signer `i` generates a secret random nonce k_i
        k_i = random.randrange(1, Q)
        signer_nonces_k[i] = k_i
        
        # Computes its public nonce share R_i = k_i * G
        R_i = (k_i * G) % Q
        signer_public_nonces_R[i] = R_i
        
        # Computes a commitment to its public nonce share: C_i = H(R_i)
        C_i = H(R_i)
        commitments_C[i] = C_i
    
    print("Each signer generates a secret nonce and broadcasts a commitment.")
    print("This round is complete.")

    # --- ROUND 2: SIGNING ---
    print("\n--- Signing Round 2: Signature Share Generation ---")
    # After receiving all commitments, signers broadcast their R_i values.
    # Each party can now verify that the received R_j matches the commitment C_j.
    # (We skip the explicit verification loop as our simulation is honest).

    # Each party computes the aggregate public nonce R and the challenge e
    aggregate_R = sum(signer_public_nonces_R.values()) % Q
    challenge_e = H(aggregate_R, public_key_P, message_m)
    print("Signers broadcast nonces, and the aggregate nonce R and challenge e are computed.")

    # Each party computes its partial signature s_i
    partial_signatures_s = {}

    def lagrange_coeff(i, S, q):
        """Computes the Lagrange coefficient for party i in set S."""
        # λ_i = product( j / (j - i) ) for j in S, j != i
        num = prod(j for j in S if j != i)
        den = prod(j - i for j in S if j != i)
        return (num * mod_inverse(den, q)) % q

    for i in signers_S:
        k_i = signer_nonces_k[i]
        x_i = private_key_shares[i]
        lambda_i = lagrange_coeff(i, signers_S, Q)
        
        # Calculate the partial signature: s_i = k_i + e * λ_i * x_i
        s_i = (k_i + challenge_e * lambda_i * x_i) % Q
        partial_signatures_s[i] = s_i

    print("Each signer has now computed their partial signature.")
    print("-" * 20)

    # --- 3. Signature Aggregation ---
    print("--- Step 3: Signature Aggregation ---")
    # A coordinator (or any participant) sums the partial signatures.
    final_s = sum(partial_signatures_s.values()) % Q
    final_signature = (aggregate_R, final_s)
    print("Partial signatures have been aggregated into the final signature (R, s).")
    print("-" * 20)

    # --- 4. Verification ---
    print("--- Step 4: Verification ---")
    print("Anyone with the group public key can verify the signature.")
    print("The verification equation is: s * G = R + e * P")
    print("where e = H(R || P || m)\n")

    # A verifier would re-calculate e using the public information
    e_verifier = H(final_signature[0], public_key_P, message_m)

    # Left-hand side of the equation
    lhs = (final_signature[1] * G) % Q
    # Right-hand side of the equation
    rhs = (final_signature[0] + e_verifier * public_key_P) % Q

    is_valid = (lhs == rhs)

    print("Final Equation Values:")
    print(f"s (final signature scalar) = {final_signature[1]}")
    print(f"G (generator)              = {G}")
    print(f"R (aggregate nonce)        = {final_signature[0]}")
    print(f"e (challenge hash)         = {e_verifier}")
    print(f"P (group public key)       = {public_key_P}")
    print("-" * 20)
    print(f"s * G (Left Side)          = {lhs}")
    print(f"R + e * P (Right Side)     = {rhs}")
    print("-" * 20)
    print(f"Signature valid: {is_valid}")

if __name__ == '__main__':
    design_threshold_signature_scheme()
<<<