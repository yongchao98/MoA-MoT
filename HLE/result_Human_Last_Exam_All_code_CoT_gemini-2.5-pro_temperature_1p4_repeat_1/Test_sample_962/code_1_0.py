import random

# A simple class to represent a participant in the FROST scheme
class Participant:
    def __init__(self, id, secret_share, t, group):
        self.id = id
        self.secret_share = secret_share
        self.group = group
        # Pre-computation for Round 1
        self.nonce_d = random.randint(1, group['q'] - 1)
        self.nonce_e = random.randint(1, group['q'] - 1)
        self.commitment_D = pow(group['g'], self.nonce_d, group['q'])
        self.commitment_E = pow(group['g'], self.nonce_e, group['q'])

def shamir_secret_sharing(t, n, secret, q):
    """Generates shares using Shamir's Secret Sharing with a trusted dealer."""
    if t > n:
        raise ValueError("Threshold t cannot be greater than the number of participants n.")
    
    # Generate a random polynomial of degree t-1
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # The secret a_0 is the group's secret key
    coeffs = [secret] + [random.randint(1, q - 1) for _ in range(t - 1)]
    
    # Generate shares for each participant i = 1, 2, ..., n
    # A share is a point on the polynomial: (i, f(i))
    shares = {}
    for i in range(1, n + 1):
        y = 0
        # Evaluate polynomial f(i)
        for j in range(t):
            y = (y + coeffs[j] * pow(i, j, q)) % q
        shares[i] = y
        
    return shares, coeffs

def lagrange_coeff(i, S, q):
    """Calculates the Lagrange coefficient L_i for a participant i in a set S."""
    num = 1
    den = 1
    for j in S:
        if i == j:
            continue
        num = (num * j) % q
        den = (den * (j - i)) % q
    # We use modular inverse for division
    return (num * pow(den, -1, q)) % q

def simple_hash(*args):
    """A simplified non-cryptographic hash function for demonstration."""
    s = "".join(map(str, args))
    return int(hash(s))

def design_frost_scheme():
    """
    Simulates the FROST 2-round threshold signature scheme.
    """
    # 1. System Parameters Setup
    # In a real system, these would be parameters of an elliptic curve group
    print("--- 1. System Setup ---")
    q = 101  # A prime order for our finite field
    g = 3    # A generator for the group
    n = 5    # Total number of participants
    t = 3    # Threshold required to sign
    group_params = {'q': q, 'g': g}
    print(f"Group Order (q): {q}, Generator (g): {g}")
    print(f"Total Participants (n): {n}, Threshold (t): {t}\n")

    # 2. Key Generation (Performed by a Trusted Dealer)
    print("--- 2. Key Generation (Trusted Dealer) ---")
    group_secret_key = random.randint(1, q - 1)
    shares, poly_coeffs = shamir_secret_sharing(t, n, group_secret_key, q)
    group_public_key = pow(g, group_secret_key, q)
    print(f"Group's master secret key (never revealed): {group_secret_key}")
    print(f"Group's public key (Y = g^secret): {group_public_key}")
    for participant_id, share in shares.items():
        print(f"  - Participant {participant_id} receives secret share: {share}")
    print("\n")
    
    # 3. Signing Protocol
    print("--- 3. Signing Protocol ---")
    message = 42 # The message to be signed
    
    # A random set of t participants will sign
    signer_ids = sorted(random.sample(list(shares.keys()), t))
    print(f"Message to sign: {message}")
    print(f"Signing participants (P): {signer_ids}\n")

    # Instantiate the participant objects for the signers
    signers = {i: Participant(i, shares[i], t, group_params) for i in signer_ids}

    # ROUND 1: Commitment
    print("--- Signing Round 1: Commitment ---")
    commitments = {}
    for pid, p_obj in signers.items():
        commitments[pid] = (p_obj.commitment_D, p_obj.commitment_E)
        print(f"Participant {pid} generates nonces and broadcasts commitments (D_{pid}, E_{pid}): ({p_obj.commitment_D}, {p_obj.commitment_E})")
    print("-> All participants broadcast their commitments to each other.\n")

    # ROUND 2: Signature Share Generation
    print("--- Signing Round 2: Signature Share Generation ---")
    # This part is executed by each participant after receiving all commitments from Round 1
    
    # First, the group commitment R is computed by all participants
    R = 1
    # Create a list of commitment tuples for hashing
    commitment_list = [cid for cid in commitments.values()]
    for j in signer_ids:
        D_j, E_j = commitments[j]
        # Binding factor rho_j = H(j, m, {D_k, E_k})
        rho_j = simple_hash(j, message, str(commitment_list)) % q
        R = (R * D_j * pow(E_j, rho_j, q)) % q
    print(f"All participants compute the same Group Commitment (R): {R}")

    # Then, the challenge c is computed by all participants
    challenge = simple_hash(group_public_key, R, message) % q
    print(f"All participants compute the same Challenge (c = H(Y, R, m)): {challenge}\n")

    # Each participant now computes their signature share
    partial_signatures = []
    print("Each participant computes their partial signature z_i = d_i + e_i*rho_i + c*L_i*s_i")
    for i in signer_ids:
        signer = signers[i]
        
        # Calculate binding factor for this signer
        rho_i = simple_hash(i, message, str(commitment_list)) % q
        
        # Calculate Lagrange coefficient for this signer
        L_i = lagrange_coeff(i, signer_ids, q)
        
        # Calculate the partial signature z_i
        term1 = signer.nonce_d
        term2 = (signer.nonce_e * rho_i) % q
        term3 = (challenge * L_i * signer.secret_share) % q
        z_i = (term1 + term2 + term3) % q
        
        partial_signatures.append(z_i)
        print(f"  - Participant {i}: z_{i} = {term1} + {term2} + {term3} (mod {q}) = {z_i}")
    print("-> All participants send their partial signature z_i to the aggregator.\n")
    
    # 4. Aggregation
    print("--- 4. Aggregation ---")
    # The aggregator sums the partial signatures to get the final signature scalar 's'
    final_s = sum(partial_signatures) % q
    print(f"Aggregator computes final signature s = sum(z_i): {final_s}")
    final_signature = (R, final_s)
    print(f"Final Threshold Signature (R, s): {final_signature}\n")

    # 5. Verification
    print("--- 5. Verification ---")
    print("Anyone can verify the signature with the group public key.")
    print("Verification equation: g^s == R * Y^c")
    
    lhs = pow(g, final_s, q)
    rhs = (R * pow(group_public_key, challenge, q)) % q
    
    # Printing each number in the final equation
    print("\n--- Final Equation Breakdown ---")
    print(f"  g = {g}")
    print(f"  s = {final_s}")
    print(f"  R = {R}")
    print(f"  Y (Group Public Key) = {group_public_key}")
    print(f"  c (Challenge) = {challenge}")
    print(f"  q = {q}\n")

    print(f"Calculating Left Hand Side (g^s mod q):")
    print(f"  {g}^{final_s} mod {q} = {lhs}")
    
    print(f"\nCalculating Right Hand Side (R * Y^c mod q):")
    print(f"  ({R} * {group_public_key}^{challenge}) mod {q} = {rhs}")

    print("\n--- Result ---")
    if lhs == rhs:
        print("Signature is VALID.")
    else:
        print("Signature is INVALID.")

# Run the simulation
if __name__ == "__main__":
    design_frost_scheme()
