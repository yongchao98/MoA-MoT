import random
import hashlib

def mod_inverse(a, m):
    """Computes the modular multiplicative inverse of a modulo m using the extended Euclidean algorithm."""
    m0, x0, x1 = m, 0, 1
    if m == 1:
        return 0
    while a > 1:
        q = a // m
        m, a = a % m, m
        x0, x1 = x1 - q * x0, x0
    if x1 < 0:
        x1 += m0
    return x1

def main():
    # 1. System Parameters
    # Set n, the total number of parties in the scheme.
    n = 5
    # Set t, the minimum number of parties required to sign (threshold).
    t = 3
    # A large prime number for our finite field operations. In a real system,
    # this would be the order of the elliptic curve group.
    prime = 1019  # A reasonably large prime for demonstration
    # The order of the group for exponent calculations (phi(prime))
    order = prime - 1

    print(f"--- System Setup ---")
    print(f"Total parties (n): {n}")
    print(f"Threshold (t): {t}")
    print(f"Finite field prime (p): {prime}")
    print(f"Group order (q): {order}\n")


    # 2. Key Generation (Simulated by a Trusted Dealer)
    # The dealer creates a secret polynomial of degree t-1.
    # The secret master key is the constant term, f(0).
    # a[0] is the master key, a[1],...,a[t-1] are other random coefficients.
    coefficients = [random.randint(1, order - 1) for _ in range(t)]
    master_secret_key = coefficients[0]

    def polynomial(x):
        res = 0
        for i in range(t):
            res = (res + coefficients[i] * (x ** i)) % order
        return res

    # Generate n secret shares for n parties. share_i = f(i)
    # Parties are indexed 1, 2, ..., n
    secret_shares = {i: polynomial(i) for i in range(1, n + 1)}

    print("--- Key Generation (Simulated DKG) ---")
    print(f"Secret polynomial created with degree {t-1}.")
    print(f"Master Secret Key (secret polynomial's constant term): {master_secret_key}")
    print("Secret shares distributed to parties:")
    for party_id, share in secret_shares.items():
        print(f"  - Party {party_id}: {share}")
    print("\n")


    # 3. Signing Protocol
    message = "Hello, world! This is a message to be signed."
    # For BLS, we hash the message to a point on the curve. Here, we hash it to an integer.
    # We use a generator 'g' for our group, let's pick a small one.
    g = 3
    # In a real system H(m) would be a complex hash-to-curve function.
    hashed_message = int(hashlib.sha256(message.encode()).hexdigest(), 16) % prime

    # Round 1: Coordinator requests signatures from a set of parties.
    # Let's assume parties 1, 3, and 5 agree to sign.
    signing_parties_indices = random.sample(list(secret_shares.keys()), t)
    signing_parties_indices.sort()
    
    print("--- Signing Protocol ---")
    print(f"Message to sign: '{message}'")
    print(f"Generator (g): {g}")
    print(f"Hashed Message (H(m)): {hashed_message}")
    print(f"Signing parties (a subset of size t={t}): {signing_parties_indices}\n")

    # Round 2: Signing parties create and send partial signatures.
    # Partial signature_i = H(m)^share_i mod p
    partial_signatures = {}
    print("--- Partial Signatures (Round 2) ---")
    for i in signing_parties_indices:
        share = secret_shares[i]
        # Partial signature is H(m) raised to the party's secret share
        partial_sig = pow(hashed_message, share, prime)
        partial_signatures[i] = partial_sig
        print(f"Party {i} (share={share}) computes partial signature: H(m)^share_{i} = {hashed_message}^{share} mod {prime} = {partial_sig}")
    print("\n")


    # 4. Signature Aggregation
    # The coordinator has the partial signatures and the indices of the signers.
    # It now computes the Lagrange coefficients to combine them.
    lagrange_coeffs = {}
    print("--- Signature Aggregation ---")
    print("Calculating Lagrange coefficients for combining signatures...")
    for i in signing_parties_indices:
        numerator = 1
        denominator = 1
        for j in signing_parties_indices:
            if i != j:
                numerator = (numerator * (0 - j)) % order
                denominator = (denominator * (i - j)) % order
        
        # We compute the final coefficient modulo the *order* of the group
        lagrange_coeffs[i] = (numerator * mod_inverse(denominator, order)) % order
        print(f"  - Lagrange coefficient L_{i}(0): {lagrange_coeffs[i]}")

    # Combine the partial signatures using the Lagrange coefficients.
    # Final signature = PRODUCT [ partial_sig_i ^ lagrange_coeff_i ] mod p
    aggregated_signature = 1
    for i in signing_parties_indices:
        term = pow(partial_signatures[i], lagrange_coeffs[i], prime)
        aggregated_signature = (aggregated_signature * term) % prime

    # Display the final equation with all the numbers
    equation_parts = []
    for i in signing_parties_indices:
        equation_parts.append(f"{partial_signatures[i]}^{lagrange_coeffs[i]}")
    
    print("\n--- Final Aggregated Signature ---")
    final_equation = " * ".join(equation_parts)
    print("The final signature is the product of each partial signature raised to its Lagrange coefficient.")
    print(f"Final Equation: Aggregated Signature = ({final_equation}) mod {prime}")
    print(f"Result: {aggregated_signature}\n")


    # 5. Verification
    # To verify, a third party checks if e(signature, g) == e(H(m), master_pk)
    # In our simulation, we check if signature == g^(H(m)*master_sk)
    # Or, equivalently, if signature == (H(m))^master_sk
    
    # We calculate what the signature *should* be if created with the master key
    correct_signature = pow(hashed_message, master_secret_key, prime)
    
    print("--- Verification ---")
    print(f"Correct signature using master key: H(m)^master_sk = {hashed_message}^{master_secret_key} mod {prime} = {correct_signature}")
    print(f"Aggregated signature from t={t} parties: {aggregated_signature}")

    if aggregated_signature == correct_signature:
        print("SUCCESS: The aggregated signature is valid.")
    else:
        print("FAILURE: The aggregated signature is invalid.")

if __name__ == '__main__':
    main()