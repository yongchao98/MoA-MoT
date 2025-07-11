import hashlib
import random

# WARNING: This is a simplified educational demonstration.
# DO NOT USE THIS CODE IN PRODUCTION.
# Real-world cryptography requires large primes, secure random number
# generation, and elliptic curve groups.

# --- 1. System Parameters (Publicly Known) ---
# In a real system, p would be a very large prime and g a generator
# of a large subgroup of Z_p*. q is the order of the subgroup.
p = 269  # A small prime number for demonstration
q = 134  # A prime order for the subgroup (p-1)/2
g = 4    # A generator of the subgroup of order q

def H( *args ):
    """A simple hash function for the scheme."""
    s = "".join(map(str, args))
    # Return an integer representation of the hash
    return int(hashlib.sha256(s.encode()).hexdigest(), 16)

def generate_lagrange_coeffs(participants, t):
    """
    Generates Lagrange coefficients for a given set of participants.
    These are needed to combine the secret shares.
    """
    coeffs = {}
    for i in participants:
        numerator = 1
        denominator = 1
        for j in participants:
            if i != j:
                numerator = (numerator * j) % q
                denominator = (denominator * (j - i)) % q
        # Compute modular inverse for the denominator
        inv_denom = pow(denominator, -1, q)
        coeffs[i] = (numerator * inv_denom) % q
    return coeffs


# --- 2. Key Generation (Performed by a Trusted Dealer) ---
def trusted_dealer_setup(n, t):
    """
    Simulates a trusted dealer generating keys for n participants with a threshold of t.
    """
    print("--- 🔑 Key Generation Phase (Trusted Dealer) ---")
    print(f"Total Participants (n): {n}, Threshold (t): {t}")

    # Generate coefficients for a polynomial of degree t-1
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    coeffs = [random.randint(1, q - 1) for _ in range(t)]
    master_secret_key = coeffs[0]  # The secret is f(0) = a_0

    # The single group public key
    group_public_key = pow(g, master_secret_key, p)

    # Generate secret shares and individual public keys for each participant
    secret_shares = {}
    public_verification_shares = {}
    for i in range(1, n + 1):
        # s_i = f(i)
        share = sum(coeffs[j] * pow(i, j, q) for j in range(t)) % q
        secret_shares[i] = share
        public_verification_shares[i] = pow(g, share, p)

    print(f"Master Secret (never shared): {master_secret_key}")
    print(f"Group Public Key (Y = g^s mod p): {group_public_key}\n")
    print("Distributed Secret Shares:")
    for i, share in secret_shares.items():
        print(f"  - Participant {i}: s_{i} = {share}")

    return master_secret_key, group_public_key, secret_shares, public_verification_shares

# --- 3. Signing Protocol ---
# This part consists of two rounds

def signing_round_one(signing_participants):
    """
    ROUND 1: Commitment. Participants generate nonces and broadcast commitments.
    This can be done before the message is known.
    """
    print("\n--- ✍️ Signing Round 1: Commitment ---")
    print(f"Participants signing: {signing_participants}")
    commitments = {}
    nonces = {}
    for i in signing_participants:
        # Each participant generates two secret nonces
        d_i, e_i = random.randint(1, q-1), random.randint(1, q-1)
        nonces[i] = (d_i, e_i)
        
        # And computes their public commitments
        D_i = pow(g, d_i, p)
        E_i = pow(g, e_i, p)
        commitments[i] = (D_i, E_i)
        print(f"  - Participant {i} generates nonces (d_{i}, e_{i}) and broadcasts commitments (D_{i}, E_{i}) = ({D_i}, {E_i})")

    return nonces, commitments

def signing_round_two(signing_participants, nonces, commitments, group_public_key, secret_shares, message):
    """
    ROUND 2: Signature Share Generation. Once message is known, participants
    compute and broadcast their partial signatures.
    """
    print("\n--- ✍️ Signing Round 2: Signature Share Generation ---")
    print(f"Message to be signed: '{message}'")
    
    # Coordinator (or each participant) forms the set of all commitments
    all_commitments_flat = [item for sublist in commitments.values() for item in sublist]
    
    # Calculate group commitment R
    R_prod = 1
    # Each participant computes the binding value and its contribution to R
    binding_values = {}
    for i in signing_participants:
        binding_input = [i, message] + all_commitments_flat
        rho_i = H(*binding_input) % q
        binding_values[i] = rho_i
        D_i, E_i = commitments[i]
        R_prod = (R_prod * D_i * pow(E_i, rho_i, p)) % p
    
    R = R_prod
    print(f"Group Commitment (R): {R}")

    # Calculate challenge c
    c = H(R, group_public_key, message) % q
    print(f"Challenge (c = H(R, Y, m)): {c}")

    # Each participant computes their partial signature
    lagrange_coeffs = generate_lagrange_coeffs(signing_participants, len(signing_participants))
    print("Lagrange Coefficients (λ): ", lagrange_coeffs)

    partial_signatures = {}
    print("Calculating partial signatures (z_i = d_i + e_i*ρ_i + c*λ_i*s_i):")
    for i in signing_participants:
        d_i, e_i = nonces[i]
        s_i = secret_shares[i]
        lambda_i = lagrange_coeffs[i]
        rho_i = binding_values[i]
        
        # This is the core FROST partial signature equation
        z_i = (d_i + e_i * rho_i + c * lambda_i * s_i) % q
        partial_signatures[i] = z_i
        print(f"  - Partial signature from participant {i} (z_{i}): {z_i}")
        
    return partial_signatures, R, c

def aggregate_signatures(partial_signatures):
    """
    The coordinator aggregates the partial signatures into one final signature value.
    """
    z = sum(partial_signatures.values()) % q
    print(f"\nAggregated signature (z = Σ z_i mod q): {z}")
    return z

# --- 4. Verification ---
def verify(group_public_key, R, z, message, c):
    """
    Anyone with the group public key can verify the signature.
    """
    print("\n--- ✅ Verification Phase ---")
    print("Verifying if g^z == R * Y^c")

    # Left side of the equation: g^z
    lhs = pow(g, z, p)

    # Right side of the equation: R * Y^c
    rhs = (R * pow(group_public_key, c, p)) % p

    print(f"Final equation components:")
    print(f"  g = {g}")
    print(f"  z = {z}")
    print(f"  R = {R}")
    print(f"  Y = {group_public_key}")
    print(f"  c = {c}")
    print(f"  p = {p}")

    print(f"\nResult:")
    print(f"  Left side (g^z mod p):  {lhs}")
    print(f"  Right side (R*Y^c mod p): {rhs}")

    if lhs == rhs:
        print("\nSignature is VALID! 🎉")
    else:
        print("\nSignature is INVALID! ☠️")

    return lhs == rhs

# --- Main Execution ---
if __name__ == "__main__":
    N_PARTICIPANTS = 3
    THRESHOLD = 2
    
    # 1. Setup
    master_sk, group_pk, shares, public_shares = trusted_dealer_setup(N_PARTICIPANTS, THRESHOLD)

    # Define which participants will sign the message (must be at least THRESHOLD)
    SIGNING_GROUP = [1, 3] # A 2-person subset of the 3 participants
    MESSAGE_TO_SIGN = "hello world"

    # 2. Signing Protocol
    # Round 1
    nonces, commitments = signing_round_one(SIGNING_GROUP)
    
    # Round 2
    partial_sigs, final_R, final_c = signing_round_two(SIGNING_GROUP, nonces, commitments, group_pk, shares, MESSAGE_TO_SIGN)
    
    # 3. Aggregation
    final_z = aggregate_signatures(partial_sigs)
    
    # The final signature is the pair (R, z)
    final_signature = (final_R, final_z)
    print(f"Final Signature (R, z): {final_signature}")

    # 4. Verification
    verify(group_pk, final_R, final_z, MESSAGE_TO_SIGN, final_c)
