import hashlib
import random
from py_ecc.secp256k1 import G, add, multiply, N

# ---
# Plan:
# 1. Set up parameters for a t-out-of-n scheme (t=3, n=5).
# 2. Simulate a trusted dealer for Distributed Key Generation (DKG).
#    - Create a secret polynomial.
#    - Distribute secret shares to n participants.
#    - Compute and publish public keys.
# 3. Perform the two-round FROST signing protocol.
#    - Round 1: Participants generate nonces and exchange commitments.
#    - Round 2: Participants calculate partial signatures and exchange them.
# 4. Aggregate the partial signatures into a final signature.
# 5. Verify the final signature and print the final verification equation.
# ---

# Step 1: Parameters and Helper Functions
N_PARTICIPANTS = 5
THRESHOLD = 3
MESSAGE = b"This is a message to be signed by the threshold group."

def sha256_as_int(data):
    """Hashes data using SHA256 and returns it as an integer modulo N."""
    return int.from_bytes(hashlib.sha256(data).digest(), 'big') % N

def point_to_bytes(point):
    """Serializes an elliptic curve point to bytes."""
    # Note: A real implementation would handle point compression and infinity.
    return point[0].to_bytes(32, 'big') + point[1].to_bytes(32, 'big')

# Step 2: Simulated Distributed Key Generation (DKG)
print("--- 1. Key Generation (Simulated with a Trusted Dealer) ---")

# The dealer creates a secret polynomial of degree t-1
# f(x) = a_0 + a_1*x + a_2*x^2 + ...
# The master secret key is a_0
coefficients = [random.randint(1, N - 1) for _ in range(THRESHOLD)]
master_secret_key = coefficients[0]
print(f"Threshold t = {THRESHOLD}, Total Participants n = {N_PARTICIPANTS}")
print(f"Dealer created a secret polynomial of degree {THRESHOLD - 1}.")

def evaluate_polynomial(x):
    """Evaluates the polynomial f(x) for a given x."""
    y = 0
    # Evaluate using Horner's method for efficiency
    for coeff in reversed(coefficients):
        y = (y * x + coeff) % N
    return y

# Generate secret shares for each participant
secret_shares = [evaluate_polynomial(i) for i in range(1, N_PARTICIPANTS + 1)]
print(f"\nSecret shares distributed to {N_PARTICIPANTS} participants.")

# Compute public keys
# Public key for the group
group_public_key = multiply(G, master_secret_key)
# Public keys for each participant's share (for verification purposes)
public_share_keys = [multiply(G, s) for s in secret_shares]
print("Group Public Key (Y) and individual public keys computed.")
print(f"Group Public Key Y = {group_public_key[0]:.4f}..., {group_public_key[1]:.4f}...")

# ---
# Step 3: Signing Protocol
# ---
print("\n--- 2. Signing Protocol ---")

# Let's assume participants 1, 3, and 5 (a threshold of 3) will sign.
signer_indices = [1, 3, 5]
signer_shares = {i: secret_shares[i-1] for i in signer_indices}
print(f"A signing group of {THRESHOLD} participants is chosen: {signer_indices}")

# --- Round 1: Commitment ---
print("\n--- Signing Round 1: Commitment ---")
nonces = {}
commitments = {}
for i in signer_indices:
    # Each signer generates a secret nonce
    nonce = random.randint(1, N - 1)
    nonces[i] = nonce
    # Each signer computes and "broadcasts" a commitment
    commitment = multiply(G, nonce)
    commitments[i] = commitment
    print(f"Participant {i}: broadcasts commitment C_{i} = {commitment[0]:.4f}..., {commitment[1]:.4f}...")

# --- Round 2: Response ---
print("\n--- Signing Round 2: Response ---")

# Each signer now has all commitments. First, they compute the group commitment R.
# Lagrange coefficient lambda_i = product(j / (j - i)) for all j in signing set
def get_lagrange_coeff(i, S):
    numerator = 1
    denominator = 1
    for j in S:
        if i == j:
            continue
        numerator = (numerator * j) % N
        denominator = (denominator * (j - i + N)) % N # Add N to avoid negative numbers
    return (numerator * pow(denominator, -1, N)) % N

# All signers compute the same group commitment R
group_commitment_R = (0, 0) # Point at infinity
for i in signer_indices:
    group_commitment_R = add(group_commitment_R, commitments[i])

print(f"All signers compute the same group commitment R = {group_commitment_R[0]:.4f}..., {group_commitment_R[1]:.4f}...")

# All signers compute the same challenge hash c
challenge_hash_c = sha256_as_int(point_to_bytes(group_commitment_R) + point_to_bytes(group_public_key) + MESSAGE)
print(f"All signers compute the challenge c = hash(R, Y, m) = {challenge_hash_c}")

# Each signer computes their partial signature z_i
partial_signatures = {}
for i in signer_indices:
    # z_i = nonce_i + (lambda_i * share_i * c)
    lambda_i = get_lagrange_coeff(i, signer_indices)
    share_contribution = (lambda_i * signer_shares[i] * challenge_hash_c) % N
    z_i = (nonces[i] + share_contribution) % N
    partial_signatures[i] = z_i
    print(f"Participant {i}: computes partial signature z_{i} = {z_i}")

# Step 4: Signature Aggregation
print("\n--- 3. Signature Aggregation ---")
# The aggregator (can be anyone) sums the partial signatures
final_signature_z = sum(partial_signatures.values()) % N
final_signature = (group_commitment_R, final_signature_z)
print(f"Aggregator computes final signature z by summing partials: z = {final_signature_z}")
print(f"Final Threshold Signature is (R, z).")

# Step 5: Signature Verification
print("\n--- 4. Signature Verification ---")

z = final_signature[1]
R = final_signature[0]
Y = group_public_key
c = challenge_hash_c # Re-compute challenge, same as before

# Verification equation: g^z == R * Y^c
# In elliptic curve notation: z*G == R + c*Y
left_side = multiply(G, z)
right_side = add(R, multiply(Y, c))

print(f"\nFinal Verification Equation: z*G = R + c*Y")
print(f"Where:")
print(f"  z = {z}")
print(f"  G is the generator point of the curve")
print(f"  R = ({R[0]}, {R[1]})")
print(f"  c = {c}")
print(f"  Y = ({Y[0]}, {Y[1]})")

print("\nExecuting verification check...")
print(f"  Left side (z*G):  ({left_side[0]}, {left_side[1]})")
print(f"  Right side (R + c*Y): ({right_side[0]}, {right_side[1]})")

# Check if the signature is valid
if left_side == right_side:
    print("\n✅ Success! The threshold signature is valid.")
else:
    print("\n❌ Failure! The threshold signature is invalid.")