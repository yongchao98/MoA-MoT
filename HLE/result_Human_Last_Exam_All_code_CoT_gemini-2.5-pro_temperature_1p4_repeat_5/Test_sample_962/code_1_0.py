import random
import hashlib

# Step 1: Cryptographic Setup (using simple integer modular arithmetic)
# In a real system, these would be parameters of an elliptic curve group like secp256k1
# Prime modulus
q = 269513333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
# Generator
g = 3

# Scheme parameters
N = 5  # Total number of parties
T = 3  # Threshold of parties required to sign

def H(data):
    """A simple hash function."""
    return int.from_bytes(hashlib.sha256(data).digest(), 'big')

def lagrange_coeff(i, S, q):
    """Calculates the Lagrange coefficient lambda_i for participant i."""
    num = 1
    den = 1
    for j in S:
        if i == j:
            continue
        num = (num * j) % q
        den = (den * (j - i + q)) % q  # Add q to ensure positive result before modulo
    return (num * pow(den, -1, q)) % q

# Step 2: Simulate DKG (Trusted Dealer Model)
print("--- 1. Key Generation (Simulated DKG) ---")
# Master secret key (unknown to any single party)
master_secret_key = random.randint(1, q - 1)
print(f"Master Secret Key (for simulation only): {master_secret_key}")

# Generate polynomial for Shamir's Secret Sharing: P(z) = a_0 + a_1*z + ... + a_{t-1}*z^{t-1}
coeffs = [master_secret_key] + [random.randint(1, q - 1) for _ in range(T - 1)]

def P(z, coeffs, q):
    """Evaluates the polynomial at a point z."""
    y = 0
    for i in range(len(coeffs)):
        y = (y + coeffs[i] * pow(z, i, q)) % q
    return y

# Generate secret shares for N parties
secret_shares = {i: P(i, coeffs, q) for i in range(1, N + 1)}
public_key_shares = {i: pow(g, s, q) for i, s in secret_shares.items()}

# The single group public key
group_public_key = pow(g, master_secret_key, q)

print(f"Group Public Key (Y): {group_public_key}\n")
print(f"Distributed {N} secret shares. Each party i holds a secret share x_i and a public share y_i = g^x_i.\n")

# --- Signing Protocol ---
print("--- 2. Signing Protocol (2 Rounds) ---")
message_to_sign = "This is a test message for a t-out-of-n signature.".encode()

# Let's assume a subset of T participants agree to sign
participants_indices = sorted(random.sample(range(1, N + 1), T))
print(f"Signing Participants (size T={T}): {participants_indices}\n")

# --- Round 1: Commitment ---
print("--- Round 1: Commitments ---")
nonces = {} # (k_i) - secret nonces
commitments = {} # (R_i) - public commitments
for i in participants_indices:
    # Each participant generates a secret nonce k_i
    k_i = random.randint(1, q - 1)
    nonces[i] = k_i
    # and computes a public commitment R_i = g^k_i
    R_i = pow(g, k_i, q)
    commitments[i] = R_i
    print(f"Participant {i}: broadcasts commitment R_{i} = {R_i}")

# --- Round 2: Signing ---
print("\n--- Round 2: Signing ---")
# Aggregate the group commitment R from all individual commitments
group_commitment_R = 1
for R_i in commitments.values():
    group_commitment_R = (group_commitment_R * R_i) % q
print(f"Aggregated Commitment (R): {group_commitment_R}")

# Compute the challenge hash e = H(R || Y || m)
challenge_e = H(group_commitment_R.to_bytes(64, 'big') + group_public_key.to_bytes(64, 'big') + message_to_sign) % q
print(f"Challenge (e): {challenge_e}")

# Each participant now computes their partial signature s_i
partial_signatures = {}
for i in participants_indices:
    # Fetch values for participant i
    k_i = nonces[i]
    x_i = secret_shares[i]
    # Calculate Lagrange coefficient for this participant
    lambda_i = lagrange_coeff(i, participants_indices, q)
    
    # Calculate partial signature: s_i = k_i + e * lambda_i * x_i
    s_i = (k_i + challenge_e * lambda_i * x_i) % q
    partial_signatures[i] = s_i
    print(f"Participant {i}: computes partial signature s_{i} = {s_i}")

# --- Signature Aggregation ---
print("\n--- 3. Signature Aggregation ---")
# The aggregator sums the partial signatures to get the final signature value s
final_s = sum(partial_signatures.values()) % q

print(f"Final Aggregated Signature is the pair (R, s):")
print(f"  R = {group_commitment_R}")
print(f"  s = {final_s}")

# --- 4. Signature Verification ---
print("\n--- 4. Verification ---")
print("Verification equation: g^s == R * Y^e")

# Left side of equation: g^s
lhs = pow(g, final_s, q)
# Right side of equation: R * Y^e
rhs = (group_commitment_R * pow(group_public_key, challenge_e, q)) % q

print(f"\nNumbers in the final equation:")
print(f"  Base (g): {g}")
print(f"  Signature scalar (s): {final_s}")
print(f"  Commitment (R): {group_commitment_R}")
print(f"  Group Public Key (Y): {group_public_key}")
print(f"  Challenge (e): {challenge_e}")

print(f"\nCalculating left side (g^s mod q):")
print(f"  {g}^{final_s} mod {q} = {lhs}")
print(f"Calculating right side (R * Y^e mod q):")
print(f"  ({group_commitment_R} * {group_public_key}^{challenge_e}) mod {q} = {rhs}")

print("\n--- Final Result ---")
if lhs == rhs:
    print("SUCCESS: Signature is valid.")
else:
    print("FAILURE: Signature is invalid.")

<<<Done>>>