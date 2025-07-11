import random
import hashlib

def design_threshold_signature_scheme():
    """
    Designs and demonstrates a 2-round t-out-of-n threshold signature scheme.
    """
    # I. System Parameters and Setup (Trusted Dealer)
    # We simulate the elliptic curve group with modular arithmetic.
    # Let G be a generator, we can treat G=1 for simplicity.
    # Scalar multiplication s * G becomes s mod p.
    # Adding points P1 + P2 becomes (p1 + p2) mod p.

    # 1. Parameters
    n = 5  # Total number of participants
    t = 3  # Threshold required to sign
    # Using a small prime for demonstration; in reality, this would be a huge 256-bit prime.
    p = 101

    # 2. Trusted Dealer generates a secret polynomial of degree t-1
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # The coefficients are chosen randomly from the finite field Z_p
    coefficients = [random.randint(1, p - 1) for _ in range(t)]

    def polynomial(x, coeffs):
        """Evaluates the polynomial at x in the finite field Z_p."""
        y = 0
        # Compute sum(coeffs[i] * x^i) mod p
        for i in range(len(coeffs)):
            y = (y + coeffs[i] * pow(x, i, p)) % p
        return y

    # 3. The master secret key is f(0)
    master_secret_key = polynomial(0, coefficients) # sk = a_0
    # The master public key is sk * G. In our simulation, PK = sk
    master_public_key = master_secret_key

    print("--- System Setup ---")
    print(f"Total participants (n): {n}")
    print(f"Threshold (t): {t}")
    print(f"Finite field prime (p): {p}")
    print(f"Secret polynomial coefficients (a_i): {coefficients}")
    print(f"Master Secret Key (sk = f(0)): {master_secret_key}")
    print(f"Master Public Key (PK): {master_public_key}\n")

    # 4. Dealer creates and distributes secret shares to n participants
    # The share for participant i is s_i = f(i)
    secret_shares = {i: polynomial(i, coefficients) for i in range(1, n + 1)}

    print("--- Secret Shares (only known to each participant) ---")
    for i, share in secret_shares.items():
        print(f"Participant {i}'s secret share (s_{i}): {share}")
    print("\n")

    # II. Two-Round Signing Protocol

    # 1. A message to be signed
    message = b"This is a very important message."

    def H(msg):
        """A hash function mapping a message to a number in Z_p."""
        # In a real BLS scheme, H would map to a point on the curve.
        # Here we map to an integer mod p.
        return int(hashlib.sha256(msg).hexdigest(), 16) % p

    message_hash = H(message)
    print("--- Signing Process ---")
    print(f"Message to sign: {message.decode()}")
    print(f"Hashed message H(m): {message_hash}\n")

    # 2. A group of t participants decides to sign
    # Let's choose participants 1, 3, and 5
    signing_participants_indices = [1, 3, 5]
    assert len(signing_participants_indices) >= t

    print(f"Participants {signing_participants_indices} will sign the message.\n")

    # 3. ROUND 1: Commitment
    # Each participant computes their partial signature and commits to it.
    # Partial signature: sigma_i = s_i * H(m)
    # Commitment: C_i = Hash(sigma_i)
    partial_signatures = {}
    commitments = {}

    print("--- Round 1: Commitments ---")
    for i in signing_participants_indices:
        s_i = secret_shares[i]
        # sigma_i = (s_i * H(m)) mod p
        sigma_i = (s_i * message_hash) % p
        partial_signatures[i] = sigma_i

        # Commit by hashing the partial signature
        comm_i = hashlib.sha256(str(sigma_i).encode()).hexdigest()
        commitments[i] = comm_i
        print(f"Participant {i} computes partial signature and broadcasts commitment: {comm_i[:16]}...")
    print("\n(The aggregator or other participants collect all commitments)\n")

    # 4. ROUND 2: Reveal & Aggregate
    # After all commitments are received, participants reveal their partial signatures.
    # The aggregator verifies the revealed signatures against the commitments.
    print("--- Round 2: Reveal & Verification ---")
    verified_partial_signatures = {}
    for i in signing_participants_indices:
        sigma_i = partial_signatures[i]
        comm_i = commitments[i]

        # Aggregator checks if the revealed signature matches the commitment
        is_valid = hashlib.sha256(str(sigma_i).encode()).hexdigest() == comm_i
        if is_valid:
            print(f"Participant {i}'s revealed signature sigma_{i} = {sigma_i} is VALID.")
            verified_partial_signatures[i] = sigma_i
        else:
            print(f"Participant {i}'s signature is INVALID. Discarding.")
    print("\n")

    # 5. Signature Reconstruction using Lagrange Interpolation
    # The final signature is sigma = sum(L_i * sigma_i) mod p
    # where L_i is the Lagrange coefficient for participant i evaluated at x=0.
    # L_i = product( j / (j - i) ) for j in signing set S, j != i
    def modular_inverse(n, prime):
        """Computes the modular inverse of n modulo prime using Fermat's Little Theorem."""
        return pow(n, prime - 2, prime)

    lagrange_coeffs = {}
    print("--- Signature Reconstruction ---")
    print("Calculating Lagrange coefficients (L_i) for x=0:")

    final_signature = 0
    sum_terms = []
    
    # Get the actual indices of participants who provided valid signatures
    valid_indices = list(verified_partial_signatures.keys())

    for i in valid_indices:
        numerator = 1
        denominator = 1
        for j in valid_indices:
            if i == j:
                continue
            # L_i = product( j / (j - i) )
            numerator = (numerator * j) % p
            denominator = (denominator * (j - i)) % p

        # L_i = numerator * mod_inverse(denominator)
        lagrange_i = (numerator * modular_inverse(denominator, p)) % p
        lagrange_coeffs[i] = lagrange_i
        print(f"  L_{i} = {lagrange_i}")

        sigma_i = verified_partial_signatures[i]
        term = (lagrange_i * sigma_i) % p
        sum_terms.append(f"({lagrange_i} * {sigma_i})")
        final_signature = (final_signature + term) % p

    equation_str = " + ".join(sum_terms)
    print(f"\nFinal Signature (sigma) = sum(L_i * sigma_i) mod {p}")
    print(f"sigma = {equation_str} mod {p}")
    print(f"sigma = {final_signature}\n")

    # III. Verification
    # A verifier checks if the signature is valid for the master public key.
    # Verification check in our simulation: sigma == PK * H(m) mod p
    print("--- Verification ---")
    expected_signature = (master_public_key * message_hash) % p
    print(f"Verifier calculates expected signature: PK * H(m) = {master_public_key} * {message_hash} mod {p} = {expected_signature}")
    print(f"Is reconstructed signature ({final_signature}) == expected signature ({expected_signature})?")

    if final_signature == expected_signature:
        print(">>> SUCCESS: The threshold signature is valid. <<<\n")
    else:
        print(">>> FAILURE: The threshold signature is invalid. <<<\n")

    # Outputting the final equation numbers as requested
    print("--- Final Equation Numbers ---")
    print(f"Reconstructed Signature: {final_signature}")
    for i in valid_indices:
        print(f"  Component from Participant {i}: L_{i} * sigma_{i} = {lagrange_coeffs[i]} * {verified_partial_signatures[i]}")
    print(f"Expected Signature: {expected_signature}")

if __name__ == '__main__':
    design_threshold_signature_scheme()