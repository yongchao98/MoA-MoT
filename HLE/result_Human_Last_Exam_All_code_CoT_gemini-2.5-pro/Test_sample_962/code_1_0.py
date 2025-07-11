import random
from typing import List, Set, Tuple

# This code requires the 'py_ecc' library.
# You can install it using: pip install py_ecc
from py_ecc.bls import G2_uncompressed as g2, hash_to_g1, sign, verify, aggregate_signatures, aggregate_pubkeys
from py_ecc.bls.ciphersuites import G2_BASIC
from py_ecc.bls.utils import privtopub
from py_ecc.typing import BLSPubkey, BLSSignature, G1Uncompressed

# Curve parameters from the BLS12-381 standard
CURVE_ORDER = 25519106679412162953284894144439160937929038198753384803433583258833938633393

# --- Helper Functions ---

def modular_inverse(a: int, m: int) -> int:
    """Computes the modular multiplicative inverse of a modulo m using Fermat's Little Theorem."""
    return pow(a, m - 2, m)

def get_lagrange_coeff(signer_id: int, signer_ids: Set[int]) -> int:
    """
    Calculates the Lagrange basis polynomial L_i(0) for a given signer_id.
    L_i(0) = product of (j / (j - i)) for all j in signer_ids where j != i.
    """
    numerator = 1
    denominator = 1
    for j in signer_ids:
        if j != signer_id:
            numerator = (numerator * j) % CURVE_ORDER
            denominator = (denominator * (j - signer_id)) % CURVE_ORDER
    return (numerator * modular_inverse(denominator, CURVE_ORDER)) % CURVE_ORDER

def aggregate_partial_signatures(
    partial_sigs: List[BLSSignature],
    lagrange_coeffs: List[int]
) -> BLSSignature:
    """
    Aggregates partial signatures using their Lagrange coefficients.
    sigma = sum(lagrange_coeff_i * partial_sig_i)
    """
    from py_ecc.bls.point_compression import decompress_g1
    from py_ecc.bls.point_manipulation import add, multiply
    from py_ecc.bls.point_compression import compress_g1

    # Decompress points to perform elliptic curve operations
    decompressed_sigs = [decompress_g1(sig) for sig in partial_sigs]
    
    # Perform scalar multiplication: lagrange_coeff * partial_sig
    scaled_sigs = [
        multiply(sig, coeff) for sig, coeff in zip(decompressed_sigs, lagrange_coeffs)
    ]

    # Sum the resulting points
    aggregated_sig_decompressed = scaled_sigs[0]
    for i in range(1, len(scaled_sigs)):
        aggregated_sig_decompressed = add(aggregated_sig_decompressed, scaled_sigs[i])

    # Compress the final point back to its signature representation
    return compress_g1(aggregated_sig_decompressed)

# --- Main Demonstration ---

def main():
    """Demonstrates a 2-round t-out-of-n threshold signature scheme."""
    
    # 1. SETUP: Define parameters for a t-out-of-n scheme
    n = 10  # Total number of participants
    t = 4   # Threshold required to sign

    print(f"--- Setting up a {t}-out-of-{n} Threshold Signature Scheme ---\n")

    # 2. SIMULATED DKG: Create keys using a trusted dealer approach
    # Create a secret polynomial of degree t-1
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # The master secret key is sk = f(0) = a_0
    poly = [random.randint(1, CURVE_ORDER - 1) for _ in range(t)]
    master_sk = poly[0]
    master_pk = privtopub(master_sk)

    # Generate secret key shares for each of the n participants
    # sk_i = f(i)
    secret_shares = {}
    participant_ids = list(range(1, n + 1))
    for i in participant_ids:
        y = 0
        for coeff_idx, coeff_val in enumerate(poly):
            y = (y + coeff_val * pow(i, coeff_idx, CURVE_ORDER)) % CURVE_ORDER
        secret_shares[i] = y

    print(f"Generated a master public key and {n} secret key shares.\n")
    
    # 3. SIGNING: A group of t participants decides to sign a message
    message = b"This is a message to be signed by the group."
    signing_participants = set(random.sample(participant_ids, t))
    
    print(f"--- Signing Protocol Start ---")
    print(f"Message to sign: '{message.decode()}'")
    print(f"Signing participants (IDs): {sorted(list(signing_participants))}\n")

    # ROUND 1: Participants create and send partial signatures
    # Each signer computes H(m) and signs it with their secret share sk_i
    message_hash_g1 = hash_to_g1(message, G2_BASIC.DST)
    
    partial_signatures = {}
    for signer_id in signing_participants:
        sk_i = secret_shares[signer_id]
        # partial_sig_i = sk_i * H(m)
        partial_sig = sign(message, sk_i, G2_BASIC)
        partial_signatures[signer_id] = partial_sig

    print("--- Round 1: Partial Signatures Created ---\n"
          "Each of the {t} participants has computed a partial signature.\n")

    # ROUND 2: Coordinator aggregates the partial signatures
    print("--- Round 2: Signature Aggregation by Coordinator ---")
    
    # For the final equation, we need the Lagrange coefficients
    coeffs_map = {
        signer_id: get_lagrange_coeff(signer_id, signing_participants)
        for signer_id in signing_participants
    }

    # Prepare data for aggregation function
    signer_list = sorted(list(signing_participants))
    sigs_to_aggregate = [partial_signatures[i] for i in signer_list]
    coeffs_to_aggregate = [coeffs_map[i] for i in signer_list]

    print("Constructing final signature: sigma = sum(L_i * partial_sig_i)\n")
    for i, signer_id in enumerate(signer_list):
        coeff = coeffs_to_aggregate[i]
        # Printing a short version of the signature for readability
        partial_sig_str = "0x" + sigs_to_aggregate[i].hex()[:10] + "..."
        print(f"Term {i+1}:")
        print(f"  Participant ID: {signer_id}")
        print(f"  Lagrange Coefficient (L_{signer_id}): {coeff}")
        print(f"  Partial Signature (sig_{signer_id}): {partial_sig_str}")
        print("-" * 20)

    # Aggregate the signatures using the computed coefficients
    final_signature = aggregate_partial_signatures(sigs_to_aggregate, coeffs_to_aggregate)

    print(f"\nFinal Aggregated Signature (sigma): 0x{final_signature.hex()}")

    # 4. VERIFICATION: Anyone can verify the final signature against the master public key
    is_valid = verify(master_pk, message, final_signature, G2_BASIC)
    
    print("\n--- Verifying Final Signature ---")
    print(f"Master Public Key used for verification.")
    print(f"Verification result: {'SUCCESS' if is_valid else 'FAILURE'}")

if __name__ == "__main__":
    main()
<<<SUCCESS>>>