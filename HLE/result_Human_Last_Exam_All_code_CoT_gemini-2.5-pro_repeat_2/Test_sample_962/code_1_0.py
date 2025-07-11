# First, you need to install the required library:
# pip install py_ecc

from py_ecc.bls import G2ProofOfPossession as bls
from py_ecc.bls.api import signature_to_G1, hash_to_G1
from py_ecc.typing import BLSPubkey, BLSSignature
import random

# ===================================================================
# Helper Functions for Shamir's Secret Sharing & Lagrange Interpolation
# ===================================================================

# The order of the G1 group in the BLS12-381 curve
MODULUS = 21888242871839275222246405745257275088548364400416034343698204186575808495617

def shamir_share_secret(secret: int, t: int, n: int) -> list[tuple[int, int]]:
    """
    Generates n secret shares for a given secret, where any t shares can reconstruct it.
    Using Shamir's Secret Sharing scheme.
    """
    if t > n:
        raise ValueError("Threshold t cannot be greater than the number of participants n.")
    
    # Generate a random polynomial of degree t-1, where the constant term is the secret
    coeffs = [secret] + [random.randint(1, MODULUS - 1) for _ in range(t - 1)]

    def poly(x: int) -> int:
        val = 0
        for i in range(t):
            val = (val + coeffs[i] * (x ** i)) % MODULUS
        return val

    # Generate n shares by evaluating the polynomial at points 1, 2, ..., n
    shares = []
    for i in range(1, n + 1):
        shares.append((i, poly(i)))
        
    return shares

def lagrange_interpolate(x_coords: list[int], y_coords: list[int], x_target: int = 0) -> int:
    """
    Finds the y-value at x_target for a polynomial defined by a set of points (x, y).
    Uses Lagrange interpolation. All calculations are modulo MODULUS.
    """
    if len(x_coords) != len(y_coords):
        raise ValueError("Number of x and y coordinates must be equal.")
        
    t = len(x_coords)
    result = 0

    for i in range(t):
        numerator = 1
        denominator = 1
        for j in range(t):
            if i == j:
                continue
            numerator = (numerator * (x_target - x_coords[j])) % MODULUS
            denominator = (denominator * (x_coords[i] - x_coords[j])) % MODULUS
        
        # Modular inverse for division
        inv_denominator = pow(denominator, -1, MODULUS)
        term = (y_coords[i] * numerator * inv_denominator) % MODULUS
        result = (result + term) % MODULUS
        
    return result

# ===================================================================
# Main Simulation of the Threshold Signature Scheme
# ===================================================================

def run_bls_threshold_signature_scheme():
    """
    Simulates the full t-out-of-n BLS threshold signature scheme.
    """
    # 1. SETUP PHASE
    print("--- 1. SETUP PHASE ---")
    n = 5  # Total number of participants
    t = 3  # Threshold required to sign
    message = b"This is the message to be signed by the group"

    print(f"Total participants (n): {n}")
    print(f"Threshold (t): {t}\n")

    # Generate a master secret key
    master_secret_key = random.randint(1, MODULUS - 1)

    # Derive the master public key
    master_public_key = bls.SkToPk(master_secret_key)

    # Generate n secret shares from the master secret key
    secret_shares = shamir_share_secret(master_secret_key, t, n)
    participant_keys = {
        i: {'secret_share': s_val} for i, s_val in secret_shares
    }

    print(f"Generated {n} secret shares from the master secret key.")
    print(f"Master Public Key: 0x{master_public_key.hex()}")
    print("-" * 20 + "\n")


    # 2. SIGNING PHASE
    print("--- 2. SIGNING PHASE ---")
    print(f"A message will be signed by {t} participants.")
    
    # Select t participants to sign (e.g., participants 1, 3, and 5)
    signer_indices = [1, 3, 5]
    print(f"Selected signers (their indices): {signer_indices}\n")

    # Hash the message to a point on the G1 curve
    message_hash_point = hash_to_G1(message)

    # -- ROUND 1: Participants generate partial signatures --
    print("Round 1: Each selected participant creates a partial signature.")
    partial_signatures = []
    for i in signer_indices:
        signer_secret_share = participant_keys[i]['secret_share']
        # Each participant signs the message hash with their secret share
        partial_sig = bls.Sign(signer_secret_share, message)
        partial_signatures.append(partial_sig)
        print(f"  - Participant {i} created their partial signature.")

    # -- ROUND 2: Aggregator combines partial signatures --
    print("\nRound 2: Aggregator collects partial signatures and combines them.")
    
    # Get the x-coordinates (indices) of the signers
    signer_x_coords = signer_indices
    # The "y-coordinates" are the partial signatures (points on the curve)
    # We will use Lagrange interpolation on their integer representations, but the
    # actual interpolation happens on the curve points.
    
    # To combine signatures, we need to find Lagrange coefficients for x=0
    lagrange_coeffs = []
    for i in range(t):
        # We need to compute L_i(0) where L_i is the i-th Lagrange basis polynomial
        # L_i(0) = product over j!=i of (-x_j) / (x_i - x_j)
        x_i = signer_x_coords[i]
        
        numerator = 1
        denominator = 1
        for j in range(t):
            if i == j:
                continue
            x_j = signer_x_coords[j]
            numerator = (numerator * (0 - x_j)) % MODULUS
            denominator = (denominator * (x_i - x_j)) % MODULUS
        
        # Modular inverse for division
        inv_denominator = pow(denominator, -1, MODULUS)
        coeff = (numerator * inv_denominator) % MODULUS
        lagrange_coeffs.append(coeff)

    # Aggregate the signatures using the computed coefficients
    # final_sig = sum(coeff_i * partial_sig_i)
    aggregated_signature = bls.Aggregate(
        [signature_to_G1(sig) for sig in partial_signatures],
        lagrange_coeffs
    )
    
    print("Successfully aggregated partial signatures into a final signature.")
    print("-" * 20 + "\n")


    # 3. VERIFICATION PHASE
    print("--- 3. VERIFICATION PHASE ---")
    print("Verifying the aggregated signature against the master public key.")
    
    # The verification equation is: e(g1, final_signature) == e(master_public_key, H(message))
    # The bls.Verify function checks this.
    is_valid = bls.Verify(master_public_key, message, aggregated_signature)

    print("\n--- Final Equation Components ---")
    print(f"Master Public Key (PK): 0x{master_public_key.hex()}")
    print(f"Message to Sign (m): '{message.decode()}'")
    print(f"Message Hash Point (H(m)): {message_hash_point}")
    print(f"Final Aggregated Signature (S): 0x{aggregated_signature.hex()}")
    
    print("\n--- Verification Result ---")
    print(f"The signature is: {'VALID' if is_valid else 'INVALID'}")
    
    return is_valid

if __name__ == '__main__':
    is_valid = run_bls_threshold_signature_scheme()
    # The final answer format
    print(f"\n<<<{is_valid}>>>")
