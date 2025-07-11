import random
import hashlib
from math import prod

# Using fastecdsa for elliptic curve operations. Install with: pip install fastecdsa
from fastecdsa.curve import secp256k1 as curve
from fastecdsa.point import Point

# --- 1. Setup Phase ---
# Let's define the parameters for our t-out-of-n scheme
N_PARTICIPANTS = 5  # Total number of participants (n)
THRESHOLD = 3      # Minimum participants required to sign (t)

# The prime order of the curve's base point G
PRIME_ORDER = curve.q

def inverse_mod(a, n):
    """Computes the modular multiplicative inverse of a modulo n."""
    return pow(a, -1, n)

def point_to_bytes(p):
    """Converts an elliptic curve point to a fixed-size byte string."""
    return p.x.to_bytes(32, 'big') + p.y.to_bytes(32, 'big')

def custom_hash(*args):
    """A helper function to hash various inputs into an integer."""
    hasher = hashlib.sha256()
    for arg in args:
        if isinstance(arg, Point):
            hasher.update(point_to_bytes(arg))
        elif isinstance(arg, str):
            hasher.update(arg.encode('utf-8'))
        elif isinstance(arg, int):
            hasher.update(arg.to_bytes(32, 'big'))
        else:
            hasher.update(bytes(arg))
    return int.from_bytes(hasher.digest(), 'big')

# --- 2. Key Generation (Simulated DKG) ---
def generate_keys(t, n, order):
    """
    Simulates a trusted dealer for key generation.
    In a real system, a Distributed Key Generation (DKG) protocol would be used.
    """
    print(f"Key Generation (Simulating DKG for {t}-out-of-{n} scheme)...")
    
    # Generate coefficients for a random polynomial of degree t-1
    # The secret master key is the constant term, poly_coeffs[0]
    poly_coeffs = [random.randrange(1, order) for _ in range(t)]
    master_secret_key = poly_coeffs[0]

    # Generate secret shares for each of the n participants
    # s_i = f(i) where f is the polynomial
    secret_shares = []
    for i in range(1, n + 1):
        share = sum(poly_coeffs[j] * (i ** j) for j in range(t)) % order
        secret_shares.append((i, share))

    # The group's public key
    group_public_key = master_secret_key * curve.G

    print("Key generation complete.\n")
    return group_public_key, secret_shares

# --- Helper for signature reconstruction ---
def get_lagrange_coeff(participant_id, signer_ids, order):
    """
    Computes the Lagrange coefficient for a given participant.
    This is essential for reconstructing the master secret from the shares.
    """
    numerator = 1
    denominator = 1
    for j in signer_ids:
        if j != participant_id:
            numerator = (numerator * j) % order
            denominator = (denominator * (j - participant_id)) % order
    
    return (numerator * inverse_mod(denominator, order)) % order


# --- 3. Signing Protocol ---
def generate_signature(group_public_key, secret_shares, t, message):
    """
    Coordinates the two-round signing protocol.
    """
    print(f"--- Starting Signature Generation for message: '{message}' ---")

    # For this example, we'll just take the first t participants as signers
    signers = secret_shares[:t]
    signer_ids = [p[0] for p in signers]
    print(f"Selected Signers (IDs): {signer_ids}\n")

    # --- ROUND 1: Commitment ---
    print("--- Round 1: Commitment ---")
    nonces = {}           # To store secret nonces k_i
    commitments = {}      # To store public commitments R_i = k_i * G

    for signer_id, _ in signers:
        k_i = random.randrange(1, PRIME_ORDER)
        nonces[signer_id] = k_i
        R_i = k_i * curve.G
        commitments[signer_id] = R_i
        print(f"  Participant {signer_id} generates a secret nonce and broadcasts commitment R_{signer_id}.")

    # Aggregate the commitments to form the group commitment R
    group_commitment_R = Point(0, 0, curve=curve) # Point at infinity
    for R_i in commitments.values():
        group_commitment_R += R_i

    print("All commitments received.\n")

    # --- ROUND 2: Signature Share ---
    print("--- Round 2: Signature Share Calculation ---")
    partial_signatures = {} # To store the signature shares z_i

    # The challenge 'c' is computed and is the same for all signers
    challenge_c = custom_hash(group_public_key, group_commitment_R, message)
    print(f"Challenge 'c' computed from group public key, group commitment, and message.")

    for signer_id, secret_share in signers:
        # Each signer calculates their partial signature
        lagrange_coeff = get_lagrange_coeff(signer_id, signer_ids, PRIME_ORDER)
        nonce_k = nonces[signer_id]
        
        # The core FROST equation for the partial signature
        z_i = (nonce_k + challenge_c * lagrange_coeff * secret_share) % PRIME_ORDER
        partial_signatures[signer_id] = z_i
        print(f"  Participant {signer_id} calculates their partial signature z_{signer_id}.")
        
    print("All partial signatures have been calculated.\n")

    # --- Aggregation ---
    print("--- Aggregating Signature ---")
    # The aggregator sums up the partial signatures
    final_signature_z = sum(partial_signatures.values()) % PRIME_ORDER
    
    # The final signature is the pair (R, z)
    final_signature = (group_commitment_R, final_signature_z)
    print("Signature aggregation complete.")
    print(f"Final Signature (R, z) successfully generated.\n")
    return final_signature


# --- 4. Verification ---
def verify_signature(group_public_key, message, signature):
    """
    Verifies the aggregated threshold signature.
    """
    R, z = signature
    print("--- Verifying Final Signature ---")

    # Re-compute the challenge 'c'
    challenge_c = custom_hash(group_public_key, R, message)

    # Verification equation: z * G == R + c * Y
    # Where Y is the group_public_key
    left_side = z * curve.G
    right_side = R + (challenge_c * group_public_key)

    print("Verification checks if: z * G == R + c * Y")
    print(f"z (final signature)   = {z}")
    print(f"G (curve generator)   = (x={curve.G.x}, y={curve.G.y})")
    print(f"R (group commitment)  = (x={R.x}, y={R.y})")
    print(f"c (challenge)         = {challenge_c}")
    print(f"Y (group public key)  = (x={group_public_key.x}, y={group_public_key.y})")
    print("\nCalculated values for the equation:")
    print(f"Left side (z * G)     = (x={left_side.x}, y={left_side.y})")
    print(f"Right side (R + c*Y)  = (x={right_side.x}, y={right_side.y})")
    
    is_valid = (left_side == right_side)
    print(f"\nSignature is valid: {is_valid}")
    return is_valid

# --- Main execution ---
if __name__ == '__main__':
    # 1. Generate keys
    group_pk, secret_shares = generate_keys(THRESHOLD, N_PARTICIPANTS, PRIME_ORDER)

    # 2. Define a message to sign
    message_to_sign = "This is a test of a two-round threshold signature."

    # 3. Generate a signature
    signature = generate_signature(group_pk, secret_shares, THRESHOLD, message_to_sign)

    # 4. Verify the signature
    verify_signature(group_pk, message_to_sign, signature)
    
    print("\n--- Testing a fraudulent message ---")
    fraudulent_message = "This is a different message."
    verify_signature(group_pk, fraudulent_message, signature)

<<<True>>>