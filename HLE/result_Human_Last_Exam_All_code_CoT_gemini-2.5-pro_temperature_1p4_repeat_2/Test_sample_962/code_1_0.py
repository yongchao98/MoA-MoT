# First, you may need to install the required library:
# pip install py_ecc

import hashlib
from typing import List, Tuple

# Import the necessary cryptographic primitives from the py_ecc library
# We use the BLS12-381 curve, a standard for BLS signatures
from py_ecc.bls12_381 import (
    G1, G2, Z1, Z2, is_inf,
    multiply,
    add,
    curve_order,
    hash_to_G2,
    pairing
)

# A simple polynomial evaluation function for f(x)
def evaluate_poly(poly: List[int], x: int) -> int:
    """Evaluates a polynomial at a given point x."""
    y = 0
    # Horner's method for polynomial evaluation
    for coeff in reversed(poly):
        y = (y * x + coeff) % curve_order
    return y

def lagrange_basis(x_coords: List[int], i_index: int, x_eval: int) -> int:
    """
    Calculates the i-th Lagrange basis polynomial at x_eval.
    x_coords are the x-coordinates of the points used for interpolation.
    i_index is the index of the basis polynomial to compute (0 to t-1).
    """
    i_x = x_coords[i_index]
    numerator = 1
    denominator = 1
    for j_index, j_x in enumerate(x_coords):
        if i_index == j_index:
            continue
        numerator = (numerator * (x_eval - j_x)) % curve_order
        denominator = (denominator * (i_x - j_x)) % curve_order

    # Modular inverse to perform division
    inv_denominator = pow(denominator, -1, curve_order)
    return (numerator * inv_denominator) % curve_order

def main():
    # 1. SETUP: Define parameters for our t-out-of-n scheme
    n = 5  # Total number of participants
    t = 3  # Threshold required to sign

    print(f"Designing a {t}-out-of-{n} BLS Threshold Signature Scheme.\n")

    # 2. KEY GENERATION (Distributed Key Generation - DKG)
    # In a real system, this is a complex protocol. Here we simulate it with a trusted dealer.
    # The dealer creates a secret polynomial of degree t-1.
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # The master secret key is f(0) = a_0
    secret_poly = [pow(i, 2, curve_order) for i in range(1, t + 1)] # Example coefficients
    master_secret_key = secret_poly[0]
    
    # The master public key is master_secret_key * G1
    master_public_key = multiply(G1, master_secret_key)

    # Each participant 'i' gets a secret share s_i = f(i)
    # We use indices 1 to n for participants
    participant_indices = list(range(1, n + 1))
    secret_shares = {i: evaluate_poly(secret_poly, i) for i in participant_indices}

    print("--- Key Generation Complete ---")
    print(f"Master secret key (never revealed): {master_secret_key}")
    print(f"Master public key (hex, compressed): {master_public_key[0].n:x}\n")

    # 3. SIGNING
    print("--- Signing Protocol (2 Rounds) ---")
    message = b"This is a very important message for the DAO to sign."
    
    # Hash the message to a point on the curve G2
    message_hash_point = hash_to_G2(message)

    # Let's assume the first `t` participants collaborate to sign
    signing_participants_indices = participant_indices[:t]
    print(f"A subset of {t} participants will sign: {signing_participants_indices}")

    # ROUND 1: Participants create their partial signatures
    # partial_sig_i = secret_share_i * H(m)
    partial_signatures = {}
    for i in signing_participants_indices:
        sk_i = secret_shares[i]
        partial_sig = multiply(message_hash_point, sk_i)
        partial_signatures[i] = partial_sig
        print(f"Participant {i} created a partial signature.")
    print("\n--- End of Signing Round 1 ---")
    
    # ROUND 2: An aggregator collects the partial signatures and combines them
    # The aggregator needs to calculate Lagrange coefficients to correctly combine the shares.
    # The final signature is sum(lambda_i * partial_sig_i) where lambda_i are Lagrange coeffs.
    
    # The aggregator calculates the coefficients at x=0 to recover the master signature
    lagrange_coeffs = {
        i: lagrange_basis(signing_participants_indices, j, 0)
        for j, i in enumerate(signing_participants_indices)
    }

    # Aggregate the signatures
    aggregated_signature = Z2 # Start with the identity element on G2
    for i in signing_participants_indices:
        lambda_i = lagrange_coeffs[i]
        partial_sig = partial_signatures[i]
        # Multiply the partial signature by its corresponding Lagrange coefficient
        term = multiply(partial_sig, lambda_i)
        # Add it to the total
        aggregated_signature = add(aggregated_signature, term)
        
    print("\n--- Start of Signing Round 2 (Aggregation) ---")
    print("Aggregator received partial signatures and is combining them...")
    print(f"Final aggregated signature (hex, compressed): {aggregated_signature[0].n:x}\n")

    # 4. VERIFICATION
    print("--- Verification ---")
    print("Anyone can verify the signature with the master public key.")
    
    # The verification equation is: e(G1, aggregated_sig) == e(master_pk, H(m))
    # where 'e' is the pairing function.
    
    # py_ecc's pairing function checks this equality.
    # pairing(pk, h(m)) * pairing(G1, -sig) == 1
    # which is equivalent to e(pk, h(m)) == e(G1, sig)
    
    p1 = pairing(master_public_key, message_hash_point)
    p2 = pairing(G1, aggregated_signature)
    
    is_valid = (p1 == p2)

    print("The final verification equation is: e(P, final_sig) == e(master_PK, H(m))")
    print("Where:")
    print(f"  P (curve generator G1):           {G1[0].n:x}")
    print(f"  final_sig (on G2):              {aggregated_signature[0].n:x}")
    print(f"  master_PK (on G1):              {master_public_key[0].n:x}")
    print(f"  H(m) (hashed message on G2):    {message_hash_point[0].n:x}\n")
    print(f"Is the signature valid? {is_valid}")
    

if __name__ == '__main__':
    main()

<<<BLS Threshold Signature Scheme>>>