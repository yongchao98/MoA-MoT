import random
from typing import List, Tuple, Dict

# py_ecc is required. You can install it with:
# pip install py_ecc
from py_ecc.bls12_381 import G1, G2, multiply, add, neg, curve_order, Z1
from py_ecc.bls import pri_to_pub, verify

def modular_inverse(a: int, m: int) -> int:
    """
    Calculates the modular inverse of a modulo m.
    Uses the Extended Euclidean Algorithm.
    """
    if m == 0:
        raise ValueError("Modulus cannot be zero")
    return pow(a, m - 2, m)

def get_lagrange_coeff(indices: List[int], i: int, field_modulus: int) -> int:
    """
    Calculates the Lagrange coefficient L_i(0) for a given set of indices.
    L_i(x) = Π_{j∈indices, j≠i} (x - j) / (i - j)
    We calculate L_i(0) = Π_{j∈indices, j≠i} (-j) / (i - j)
    """
    num = 1
    den = 1
    for j in indices:
        if i == j:
            continue
        num = (num * (-j)) % field_modulus
        den = (den * (i - j)) % field_modulus
    
    return (num * modular_inverse(den, field_modulus)) % field_modulus

def trusted_dealer_dkg(n: int, t: int) -> Tuple[int, 'G1', Dict[int, int], Dict[int, 'G1']]:
    """
    Simulates a trusted dealer for Distributed Key Generation (DKG).
    Generates a master key and n private shares for a t-out-of-n scheme.
    
    Args:
        n: Total number of participants.
        t: Threshold required for signing.
    
    Returns:
        A tuple containing:
        - master_sk (int): The master secret key.
        - master_pk (G1 point): The master public key.
        - private_shares (dict): A dictionary mapping participant index to their private key share.
        - public_shares (dict): A dictionary mapping participant index to their public verification key.
    """
    if not (1 < t <= n):
        raise ValueError("Threshold t must be between 2 and n")
        
    # Generate a random polynomial of degree t-1
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    coefficients = [random.randint(1, curve_order - 1) for _ in range(t)]
    polynomial = lambda x: sum(coefficients[j] * (x ** j) for j in range(t)) % curve_order
    
    # Master secret key is f(0)
    master_sk = polynomial(0)
    master_pk = pri_to_pub(master_sk) # master_sk * G1
    
    # Generate private and public shares for n participants
    # Participant indices are 1, 2, ..., n
    private_shares = {i: polynomial(i) for i in range(1, n + 1)}
    public_shares = {i: pri_to_pub(sk) for i, sk in private_shares.items()}
    
    return master_sk, master_pk, private_shares, public_shares

def main():
    """Main function to demonstrate the 2-round threshold signature scheme."""
    
    # --- 0. Parameters ---
    N_PARTICIPANTS = 5
    THRESHOLD = 3
    MESSAGE = b"This is a message for the two-round threshold signature scheme."
    
    print(f"Scheme Parameters: n={N_PARTICIPANTS}, t={THRESHOLD}\n")

    # --- SETUP: Key Generation (done once by a trusted dealer) ---
    print("--- SETUP: Running Trusted Dealer for Key Generation ---")
    master_sk, master_pk, private_shares, public_shares = trusted_dealer_dkg(N_PARTICIPANTS, THRESHOLD)
    print("Key generation complete. Shares distributed to participants.\n")

    # --- SIGNING PROTOCOL ---
    
    # A random subset of t participants will sign the message
    signer_indices = sorted(random.sample(list(private_shares.keys()), THRESHOLD))
    print(f"--- SIGNING: A group of {THRESHOLD} participants will sign: {signer_indices} ---")
    
    # Hash the message to a point on the curve G2
    # This is H(m) in the BLS scheme
    msg_hash_point = G2.hash_to_curve(MESSAGE)
    
    # --- Round 1: BROADCAST ---
    # Each signer computes and broadcasts their partial signature
    print("\n--- Round 1: Signers compute and broadcast partial signatures ---")
    partial_signatures = {}
    for i in signer_indices:
        sk_i = private_shares[i]
        # Partial signature: sigma_i = sk_i * H(m)
        partial_sig = multiply(msg_hash_point, sk_i)
        partial_signatures[i] = partial_sig
        print(f"Participant {i} has computed and broadcasted their partial signature.")

    # --- Round 2: AGGREGATION ---
    # A combiner receives the partial signatures and creates the final signature.
    print("\n--- Round 2: A combiner receives and aggregates signatures ---")
    
    # For robustness, the combiner should verify each partial signature.
    # verify(e(PK_i, H(m)) == e(G1, sigma_i)). We will skip this for brevity.
    
    print("Calculating Lagrange coefficients for the final equation...")
    lagrange_coeffs = {i: get_lagrange_coeff(signer_indices, i, curve_order) for i in signer_indices}
    
    # Combine the signatures using the Lagrange coefficients
    # final_sigma = SUM(lagrange_coeff_i * partial_sig_i)
    signature_components = []
    for i in signer_indices:
        coeff = lagrange_coeffs[i]
        partial_sig = partial_signatures[i]
        # component = L_i * sigma_i
        component = multiply(partial_sig, coeff)
        signature_components.append(component)
    
    # Add the components together to get the final signature point
    final_signature = signature_components[0]
    for i in range(1, len(signature_components)):
        final_signature = add(final_signature, signature_components[i])

    print("\n--- The Final Signature Equation ---")
    equation_str = "Final Signature = "
    for i, idx in enumerate(signer_indices):
        equation_str += f"({lagrange_coeffs[idx]}) * (Partial Sig from {idx})"
        if i < len(signer_indices) - 1:
            equation_str += " + "
    print(equation_str)
    print("\nSignature aggregation complete.")
    
    # --- VERIFICATION ---
    print("\n--- VERIFICATION: Verifying the final signature against the master public key ---")
    is_valid = verify(MESSAGE, master_pk, final_signature)
    
    print(f"Message: {MESSAGE.decode()}")
    print(f"Verification successful: {is_valid}")
    
    if is_valid:
        print("The 2-round threshold signature is valid!")
    else:
        print("The threshold signature is INVALID!")

if __name__ == "__main__":
    main()
<<<True>>>