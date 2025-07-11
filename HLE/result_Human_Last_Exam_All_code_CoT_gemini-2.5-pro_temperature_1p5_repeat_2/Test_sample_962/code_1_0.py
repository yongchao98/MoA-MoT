import secrets
# This code requires the py_ecc library. Please install it using:
# pip install py_ecc
from py_ecc.bls12_381 import G1, G2, curve_order, pairing, add, multiply, hash_to_G2

# Helper function for modular inverse, used for Lagrange coefficients
def mod_inverse(value, modulus):
    """Calculates the modular multiplicative inverse of value mod modulus."""
    return pow(value, -1, modulus)

def generate_polynomial(degree):
    """Generates a random polynomial of a given degree."""
    # Coefficients are chosen from the field Z_q where q is the curve order
    return [secrets.randbelow(curve_order) for _ in range(degree + 1)]

def evaluate_polynomial(poly, x):
    """Evaluates a polynomial f(x) at a given point x."""
    if x == 0:
        return poly[0]
    result = 0
    # Evaluate using Horner's method for efficiency
    for coeff in reversed(poly):
        result = (result * x + coeff) % curve_order
    return result

def lagrange_coefficient(i, signers_indices, modulus):
    """
    Calculates the Lagrange basis polynomial l_i(0) for a signer i.
    l_i(0) = product_{j in S, j!=i} (j / (j-i))
    """
    numerator = 1
    denominator = 1
    for j in signers_indices:
        if i == j:
            continue
        numerator = (numerator * j) % modulus
        denominator = (denominator * (j - i)) % modulus
    
    # We need the modular inverse of the denominator
    inv_denominator = mod_inverse(denominator, modulus)
    return (numerator * inv_denominator) % modulus

def design_bls_tss(t, n, signers_indices):
    """
    Designs and demonstrates a (t,n) BLS Threshold Signature Scheme.
    """
    print("--- 1. Scheme Setup ---")
    print(f"Designing a {t}-out-of-{n} threshold signature scheme.\n")
    if len(signers_indices) < t:
        print(f"Error: Not enough signers. Have {len(signers_indices)}, need at least {t}.")
        return

    # === Key Generation (Performed by a trusted dealer or DKG protocol) ===
    # A secret polynomial of degree t-1
    poly = generate_polynomial(t - 1)
    
    # Master secret key is f(0)
    master_secret_key = poly[0]
    
    # Master public key is g1 * msk
    master_public_key = multiply(G1, master_secret_key)
    
    print(f"Master Secret Key (f(0)): {master_secret_key} (private, for demo only)")
    print(f"Master Public Key (g*f(0)): {master_public_key}\n")
    
    # Generate secret shares and public verification keys for all n participants
    secret_shares = {}
    public_shares = {}
    print("Generating keys for all participants...")
    for i in range(1, n + 1):
        sk_i = evaluate_polynomial(poly, i)
        pk_i = multiply(G1, sk_i)
        secret_shares[i] = sk_i
        public_shares[i] = pk_i
        # print(f"  Participant {i}: sk_{i}={sk_i}, pk_{i}={pk_i}") # Optional: can be noisy
    
    print("\n--- 2. Signing Protocol (Two Rounds) ---")
    message = b"This is the message to be signed by the threshold group"
    print(f"Message to sign: '{message.decode()}'")
    
    # Hash message to a point on the G2 curve
    H_m = hash_to_G2(message)
    print(f"H(message) on G2: {H_m}\n")
    
    print(f"A group of {len(signers_indices)} participants will sign: {signers_indices}")
    print("Round 1: Signers broadcast their commitment to sign this message.")
    print("Round 2: After receiving commitments, signers compute and broadcast partial signatures.\n")
    
    # === Partial Signature Generation ===
    partial_signatures = {}
    print("Calculating Partial Signatures (sigma_i = H(m) * sk_i):")
    for i in signers_indices:
        sk_i = secret_shares[i]
        partial_sig = multiply(H_m, sk_i)
        partial_signatures[i] = partial_sig
        print(f"  Participant {i}: Partial Sig = {partial_sig}")
        
    print("\n--- 3. Signature Aggregation ---")
    
    # === Aggregation ===
    # An aggregator collects the partial signatures and combines them.
    # Calculate Lagrange coefficients for each signing participant
    lagrange_coeffs = {
        i: lagrange_coefficient(i, signers_indices, curve_order) 
        for i in signers_indices
    }
    
    # Combine signatures using the coefficients
    # final_signature = SUM( L_i * sigma_i )
    final_signature_components = []
    for i in signers_indices:
        L_i = lagrange_coeffs[i]
        term = multiply(partial_signatures[i], L_i)
        final_signature_components.append(term)
        
    final_signature = final_signature_components[0]
    for i in range(1, len(final_signature_components)):
        final_signature = add(final_signature, final_signature_components[i])
        
    print("The final signature is constructed using the formula:")
    print("  Final_Signature = SUM_{i in S} (LagrangeCoeff_i * PartialSignature_i)\n")
    print("With the following values:")
    equation_str = "Final Signature = "
    for i in signers_indices:
        print(f"  Participant {i}: Lagrange Coeff L_{i} = {lagrange_coeffs[i]}")
        equation_str += f"{lagrange_coeffs[i]} * PartialSig_{i} + "
    print(f"\nEquation: {equation_str[:-3]}") # Remove trailing " + "
    
    print(f"\nFinal Aggregated Signature on G2: {final_signature}")

    print("\n--- 4. Verification ---")
    print("To verify, we check if the pairing equation holds:")
    print("  e(G1, Final_Signature) == e(Master_Public_Key, H(m))\n")
    
    # Verification check
    p1 = pairing(G1, final_signature)
    p2 = pairing(master_public_key, H_m)
    
    print(f"Left side (e(G1, Final_Signature)):  {p1}")
    print(f"Right side (e(MPK, H(m))):        {p2}")
    
    is_valid = (p1 == p2)
    print(f"\nVerification successful: {is_valid}")
    return is_valid

if __name__ == '__main__':
    # --- Parameters ---
    # t: threshold number of signers required
    # n: total number of participants in the group
    T_THRESHOLD = 3
    N_TOTAL_PARTICIPANTS = 5
    
    # A subset of participants who will cooperate to sign the message
    # Must have at least 't' members
    signing_participants_indices = [1, 3, 5]

    result = design_bls_tss(T_THRESHOLD, N_TOTAL_PARTICIPANTS, signing_participants_indices)
    
    # The final output is the boolean result of the verification
    # Using the specified format for the final answer.
    print(f"\n<<<final_answer>>>")
    print(f"<<<{result}>>>")