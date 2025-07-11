import hashlib
import secrets
import random
from ecdsa import SECP256k1
from ecdsa.util import number_to_string, string_to_number

# --- Elliptic Curve Parameters ---
CURVE = SECP256k1
G = CURVE.generator
ORDER = CURVE.order

# --- Scheme Parameters ---
N_PARTIES = 5  # Total number of parties in the group
THRESHOLD = 3  # Minimum number of parties required to sign

def mod_inverse(k, n):
    """Computes the modular inverse of k modulo n using the extended Euclidean algorithm."""
    return pow(k, -1, n)

def polynomial_eval(coeffs, x, prime):
    """Evaluates a polynomial with coefficients `coeffs` at point `x`."""
    y = 0
    for coeff in reversed(coeffs):
        y = (y * x + coeff) % prime
    return y

def lagrange_coeff(party_id, signers, prime):
    """
    Calculates the Lagrange coefficient for a given party_id within a set of signers.
    This is used to reconstruct the secret from shares.
    The formula is L_i(0) = product of (j / (j - i)) for all j in signers where j != i.
    """
    num = 1
    den = 1
    for j in signers:
        if party_id == j:
            continue
        num = (num * j) % prime
        den = (den * (j - party_id)) % prime
    
    return (num * mod_inverse(den, prime)) % prime

def hash_to_int(*args):
    """Hashes one or more arguments to an integer modulo the curve order."""
    hasher = hashlib.sha256()
    for arg in args:
        if isinstance(arg, ecdsa.ellipticcurve.Point):
            hasher.update(arg.to_string('compressed'))
        elif isinstance(arg, int):
            hasher.update(arg.to_bytes(32, 'big'))
        elif isinstance(arg, bytes):
            hasher.update(arg)
        else:
            raise TypeError(f"Unsupported type for hashing: {type(arg)}")
    return string_to_number(hasher.digest()) % ORDER

def simulate_dkg(t, n):
    """
    Simulates a Distributed Key Generation (DKG) using a trusted dealer.
    In a real system, this would be an interactive protocol among the parties.
    """
    print(f"=== Simulating DKG for {t}-out-of-{n} scheme ===")
    # 1. Create a secret polynomial of degree t-1
    # The secret key `x` is the constant term (coeffs[0])
    poly_coeffs = [secrets.randbelow(ORDER) for _ in range(t)]
    secret_key_x = poly_coeffs[0]
    
    # 2. Generate secret shares for each of the n parties
    # Each party `i` gets a share `x_i = P(i)`
    secret_shares = {i: polynomial_eval(poly_coeffs, i, ORDER) for i in range(1, n + 1)}
    
    # 3. The group public key is Y = x * G
    group_public_key = secret_key_x * G
    
    print(f"Secret Key (x): [REDACTED]")
    print(f"Group Public Key (Y): {group_public_key.to_string('compressed').hex()}")
    print("-" * 20)
    return secret_shares, group_public_key

def main():
    """Main function to orchestrate the threshold signature scheme."""
    
    # --- 0. Setup ---
    message = b"This is a message to be signed by the threshold group."
    secret_shares, group_public_key = simulate_dkg(THRESHOLD, N_PARTIES)
    
    # Randomly choose a set of `t` signers from the `n` parties
    all_party_ids = list(range(1, N_PARTIES + 1))
    signer_ids = sorted(random.sample(all_party_ids, THRESHOLD))
    
    print(f"A group of {THRESHOLD} signers will sign the message.")
    print(f"Participant IDs: {signer_ids}")
    print(f"Message to sign: '{message.decode()}'")
    print("-" * 20)

    # --- 1. Signing: Round 1 (Commitment) ---
    print("\n=== Signing Round 1: Commitment ===")
    # Each signer `i` generates a secret nonce `k_i` and a public commitment `R_i = k_i * G`
    nonces = {i: secrets.randbelow(ORDER) for i in signer_ids}
    commitments = {i: nonces[i] * G for i in signer_ids}
    
    for i in signer_ids:
        print(f"Party {i}:")
        print(f"  - Secret Nonce (k_{i}): [REDACTED]")
        print(f"  - Public Commitment (R_{i}): {commitments[i].to_string('compressed').hex()}")
    
    # --- 2. Signing: Round 2 (Signature Share Generation) ---
    print("\n=== Signing Round 2: Signature Share Generation ===")
    # Each party receives all commitments and computes the group commitment `R`
    group_commitment_R = sum(commitments.values(), CURVE.generator * 0)
    print(f"Group Commitment (R = sum(R_i)): {group_commitment_R.to_string('compressed').hex()}")
    
    # Each party computes the same challenge `e = H(R, Y, m)`
    challenge_e = hash_to_int(group_commitment_R, group_public_key, message)
    print(f"Challenge (e = H(R, Y, m)): {challenge_e}")
    
    # Each signer `i` computes their partial signature `s_i = k_i + e * lambda_i * x_i`
    partial_signatures = {}
    for i in signer_ids:
        signer_secret_share = secret_shares[i]
        lambda_i = lagrange_coeff(i, signer_ids, ORDER)
        s_i = (nonces[i] + challenge_e * lambda_i * signer_secret_share) % ORDER
        partial_signatures[i] = s_i
        print(f"Party {i}:")
        print(f"  - Lagrange Coeff (lambda_{i}): {lambda_i}")
        print(f"  - Partial Signature (s_{i}): {s_i}")

    # --- 3. Aggregation ---
    print("\n=== Aggregation ===")
    # The aggregator sums the partial signatures to get the final signature `s`
    final_signature_s = sum(partial_signatures.values()) % ORDER
    
    print(f"Final aggregated signature (s = sum(s_i)): {final_signature_s}")
    print("\nFinal Signature is the tuple (R, s):")
    print(f"  R: {group_commitment_R.to_string('compressed').hex()}")
    print(f"  s: {final_signature_s}")
    
    # --- 4. Verification ---
    print("\n=== Verification ===")
    print("A verifier needs only the group public key (Y), the message (m), and the signature (R, s).")
    
    # The verifier computes the left and right sides of the Schnorr equation:
    # s * G == R + e * Y
    
    # Re-compute challenge `e` just as the signers did
    verifier_challenge_e = hash_to_int(group_commitment_R, group_public_key, message)
    
    left_side = final_signature_s * G
    right_side = group_commitment_R + verifier_challenge_e * group_public_key
    
    print("\nVerification Equation: s * G == R + e * Y")
    print("-" * 40)
    print(f"s = {final_signature_s}")
    print(f"G = ({G.x()}, {G.y()})")
    print(f"R = ({group_commitment_R.x()}, {group_commitment_R.y()})")
    print(f"e = {verifier_challenge_e}")
    print(f"Y = ({group_public_key.x()}, {group_public_key.y()})")
    print("-" * 40)
    print(f"Left side (s * G)  : {left_side.to_string('compressed').hex()}")
    print(f"Right side (R + e*Y): {right_side.to_string('compressed').hex()}")
    
    if left_side == right_side:
        print("\n---> SUCCESS: The signature is valid! <---")
    else:
        print("\n---> FAILURE: The signature is invalid! <---")

if __name__ == '__main__':
    main()