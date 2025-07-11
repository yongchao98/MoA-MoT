import os
import hashlib
from fastecdsa.curve import secp256k1
from fastecdsa.point import Point

# Use the secp256k1 curve, common in cryptocurrencies like Bitcoin and Ethereum
C = secp256k1
# G is the generator point of the curve
G = C.G
# q is the order of the curve's base field
q = C.q

def int_to_hex(i):
    """Helper to convert an integer to a hex string for hashing."""
    return hex(i)[2:]

def point_to_hex(p):
    """Helper to convert an elliptic curve point to a hex string for hashing."""
    return int_to_hex(p.x) + int_to_hex(p.y)

def challenge_hash(public_key, group_commitment, message):
    """Computes the challenge `c` by hashing the public key, commitment, and message."""
    h = hashlib.sha25le()
    h.update(point_to_hex(public_key).encode('utf-8'))
    h.update(point_to_hex(group_commitment).encode('utf-8'))
    h.update(message.encode('utf-8'))
    return int(h.hexdigest(), 16) % q

def shamir_share(secret, threshold, num_parties):
    """
    Generates Shamir's secret shares for a given secret.
    Returns a list of (x, y) pairs where y = f(x) for a random polynomial f.
    """
    if threshold > num_parties:
        raise ValueError("Threshold cannot be greater than the number of parties.")
    
    # Generate a random polynomial of degree (threshold - 1)
    # The constant term is the secret
    coeffs = [secret] + [int.from_bytes(os.urandom(32), 'big') % q for _ in range(threshold - 1)]

    def poly(x):
        val = 0
        # Evaluate polynomial using Horner's method for efficiency
        for coeff in reversed(coeffs):
            val = (val * x + coeff) % q
        return val

    # Generate points on the polynomial for each party
    shares = []
    for i in range(1, num_parties + 1):
        shares.append((i, poly(i)))
        
    return shares

def get_lagrange_coeffs(participant_indices):
    """
    Computes the Lagrange coefficients for a given set of participants.
    These are used to combine the signature shares correctly.
    """
    coeffs = {}
    for i in participant_indices:
        numerator = 1
        denominator = 1
        for j in participant_indices:
            if i != j:
                numerator = (numerator * j) % q
                # We need modular inverse for division
                denominator = (denominator * (j - i)) % q
        coeffs[i] = (numerator * pow(denominator, -1, q)) % q
    return coeffs

# --- Main Protocol ---
print("--- Setup ---")
t = 3  # Threshold
n = 5  # Total number of parties
message = "This is a test of a two-round threshold signature scheme"

print(f"Configuration: t={t}, n={n}\n")

# 1. KEY GENERATION (Trusted Dealer simulation)
print("1. Key Generation")
master_private_key = int.from_bytes(os.urandom(32), 'big') % q
master_public_key = G * master_private_key
print(f"Master Public Key (Y): ({master_public_key.x}, {master_public_key.y})")

# Split the master key into n shares
all_shares = shamir_share(master_private_key, t, n)
# Each party i has a private share s_i = shares[i-1][1]
# We'll just keep the list for this simulation
print(f"Generated {n} private key shares.\n")

# 2. SIGNING PROTOCOL
print("2. Signing Protocol")
# Let's assume the first t parties participate
# Participant indices are 1-based, so we take 1, 2, 3
participant_indices = [p[0] for p in all_shares[:t]]
participant_shares = {p[0]: p[1] for p in all_shares[:t]}

print(f"Signing participants (indices): {participant_indices}")
print(f"Message to sign: '{message}'\n")

# --- ROUND 1: Commitment ---
print("--- Round 1: Commitments ---")
nonces = {} # To store k_i for each participant i
commitments = {} # To store R_i for each participant i

for i in participant_indices:
    # Each participant generates a secret nonce k_i
    nonces[i] = int.from_bytes(os.urandom(32), 'big') % q
    # And computes their public commitment R_i = G * k_i
    commitments[i] = G * nonces[i]
    print(f"Participant {i}: Generated commitment R_{i}")

# Participants broadcast their commitments to each other.
# We simulate this by simply moving to the next step.
print("\n--- Round 2: Signature Shares ---")

# All participants compute the group commitment R
group_commitment_R = Point(0, 0, curve=C) # Point at infinity
for R_i in commitments.values():
    group_commitment_R += R_i
print(f"Group Commitment (R): ({group_commitment_R.x}, {group_commitment_R.y})")

# All participants compute the same challenge hash c
c = challenge_hash(master_public_key, group_commitment_R, message)
print(f"Challenge (c): {c}")

# All participants compute the Lagrange coefficients
lagrange_coeffs = get_lagrange_coeffs(participant_indices)
print(f"Lagrange Coefficients (lambda): {lagrange_coeffs}\n")

# Each participant computes their partial signature z_i
signature_shares = {}
for i in participant_indices:
    lambda_i = lagrange_coeffs[i]
    k_i = nonces[i]
    s_i = participant_shares[i]
    # z_i = k_i + c * lambda_i * s_i
    z_i = (k_i + c * lambda_i * s_i) % q
    signature_shares[i] = z_i
    print(f"Participant {i}: Generated signature share z_{i} = {z_i}")


# 3. SIGNATURE AGGREGATION
print("\n3. Signature Aggregation")
# A coordinator sums the signature shares to get the final signature scalar z
final_z = sum(signature_shares.values()) % q

print(f"Aggregated signature scalar (z): {final_z}")
print(f"Final Signature is the pair (R, z)\n")


# 4. VERIFICATION
print("4. Verification")
print("Anyone with the master public key can verify the signature.")
print("Verification equation: G * z == R + Y * c\n")

# The verifier calculates the left-hand side (LHS) of the equation
lhs = G * final_z
print(f"Left Side (G * z):")
print(f"  x = {lhs.x}")
print(f"  y = {lhs.y}\n")

# The verifier calculates the right-hand side (RHS) of the equation
rhs = group_commitment_R + master_public_key * c
print(f"Right Side (R + Y * c):")
print(f"  x = {rhs.x}")
print(f"  y = {rhs.y}\n")

# Final check
assert lhs == rhs, "Signature verification FAILED!"
print("Signature is VALID because the left and right sides of the equation match.")
print("\nFinal Equation Output:")
print(f"({lhs.x}, {lhs.y}) == ({rhs.x}, {rhs.y})")

<<<VALID>>>