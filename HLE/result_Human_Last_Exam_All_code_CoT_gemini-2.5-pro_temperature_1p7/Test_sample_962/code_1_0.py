import random
import hashlib

# --- Cryptographic Primitives (Simplified) ---

# This is a large prime number for our finite field operations.
# In a real implementation, this would be the order of an elliptic curve group.
Q = 2**256 - 2**32 - 977 

# A generator for our group. Just a number in this simulation.
G = 3 

def extended_gcd(a, b):
    """Extended Euclidean Algorithm to find modular inverse."""
    if a == 0:
        return b, 0, 1
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y

def mod_inverse(n, q):
    """Calculate the modular inverse of n modulo q."""
    _, inv, _ = extended_gcd(n, q)
    return inv % q

def get_lagrange_coeff(i, S, q):
    """
    Calculate the Lagrange coefficient for participant i in the set S.
    L_i = product of j/(j-i) for all j in S where j != i.
    """
    S_list = list(S)
    numerator = 1
    denominator = 1
    for j in S_list:
        if i == j:
            continue
        numerator = (numerator * j) % q
        denominator = (denominator * (j - i)) % q
    return (numerator * mod_inverse(denominator, q)) % q

def H(*args):
    """A simple hash function simulation."""
    hasher = hashlib.sha256()
    for arg in args:
        hasher.update(str(arg).encode())
    return int(hasher.hexdigest(), 16) % Q

# --- Simulation Setup ---

# Parameters for the t-out-of-n scheme
N_PARTIES = 5
THRESHOLD = 3

print(f"--- System Setup: A {THRESHOLD}-out-of-{N_PARTIES} Threshold Signature Scheme ---\n")

# --- 1. Key Generation (Dealer Simulation) ---
# In a real system, this is done with a Distributed Key Generation (DKG) protocol.
# No single party ever knows the master secret_key.

# Generate coefficients for a polynomial of degree t-1
# f(x) = secret_key + a1*x + a2*x^2 + ...
secret_key = random.randint(1, Q - 1)
coeffs = [secret_key] + [random.randint(1, Q - 1) for _ in range(THRESHOLD - 1)]

# Generate key shares for each of the N parties
secret_shares = {}
public_key_shares = {}
for i in range(1, N_PARTIES + 1):
    # Evaluate polynomial f(i) to get the secret share s_i
    s_i = sum(c * (i ** j) for j, c in enumerate(coeffs)) % Q
    secret_shares[i] = s_i
    # The public part of the share is P_i = s_i * G
    public_key_shares[i] = pow(G, s_i, Q)

# The master public key is Y = secret_key * G
master_public_key = pow(G, secret_key, Q)

print("1. Key Generation Complete:")
print(f"   - Master Public Key (Y): {master_public_key}")
print(f"   - {N_PARTIES} secret shares distributed. Each party `i` has a secret `s_i` and public key share `Y_i`.\n")

# --- 2. Signing Protocol ---

# Define the message to be signed
message = "This is a test message for the FROST signing protocol"
print(f"2. Signing a Message:")
print(f"   - Message (m): '{message}'")

# A subset of THRESHOLD participants will collaborate to sign
signing_participants = set(range(1, THRESHOLD + 1))
print(f"   - Signing Participants (S): {list(signing_participants)}\n")

# -- Round 1: Commitment --
print("--- Round 1: Commitment ---")
nonces = {}
commitments = {}
for i in signing_participants:
    # Each participant generates a secret nonce k_i
    k_i = random.randint(1, Q-1)
    nonces[i] = k_i
    
    # And creates a public commitment R_i = k_i * G
    R_i = pow(G, k_i, Q)
    commitments[i] = R_i
    print(f"   - Participant {i}: Generates secret nonce and broadcasts commitment R_{i} = {R_i}")

print("\n--- Round 2: Signing ---")
# Participants broadcast their nonces after seeing all commitments.
# An aggregator (can be one of the participants) computes the group commitment and challenge.

# The group commitment R is the product of individual commitments
group_commitment_R = 1
for i in signing_participants:
    group_commitment_R = (group_commitment_R * commitments[i]) % Q
print(f"   - Aggregator computes group commitment R: {group_commitment_R}")
    
# The challenge `c` is a hash of the master public key, group commitment, and message.
challenge_c = H(master_public_key, group_commitment_R, message)
print(f"   - Aggregator computes challenge c = H(Y, R, m): {challenge_c}")

# Each participant now computes their partial signature z_i
partial_signatures = {}
for i in signing_participants:
    # Lagrange coefficient L_i for this participant
    lagrange_coeff_i = get_lagrange_coeff(i, signing_participants, Q)
    
    # Partial signature: z_i = k_i + c * L_i * s_i
    s_i = secret_shares[i]
    k_i = nonces[i]
    z_i = (k_i + challenge_c * lagrange_coeff_i * s_i) % Q
    partial_signatures[i] = z_i
    print(f"   - Participant {i}: Computes and sends partial signature z_{i} = {z_i}")

# -- Aggregation --
print("\n--- 3. Signature Aggregation ---")
# The final signature is the sum of the partial signatures
final_signature_z = sum(partial_signatures.values()) % Q

print(f"   - Aggregator combines partials into the final signature z = sum(z_i): {final_signature_z}")

# The complete threshold signature is the pair (R, z)
print(f"\n=> The final Schnorr Threshold Signature is (R, z):")
print(f"   R = {group_commitment_R}")
print(f"   z = {final_signature_z}\n")


# --- 4. Signature Verification ---
print("--- 4. Verification ---")
print("A verifier only needs the message, the master public key Y, and the signature (R, z).")
print("Verification Equation: z * G == R * (Y^c)")

# Left side of equation: z * G
lhs = pow(G, final_signature_z, Q)

# Right side of equation: R * (Y^c)
rhs = (group_commitment_R * pow(master_public_key, challenge_c, Q)) % Q

print("\nFinal Equation Components:")
print(f"   z*G (LHS) : {lhs}")
print(f"   R*(Y^c) (RHS) : {rhs}")

# Check if the signature is valid
is_valid = (lhs == rhs)

print(f"\nResult: The signature is {'VALID' if is_valid else 'INVALID'}.")

# --- Why it works (Mathematical Proof) ---
# z * G = (sum(k_i + c * L_i * s_i)) * G
#       = (sum(k_i) + c * sum(L_i * s_i)) * G
#       = sum(k_i) * G + c * sum(L_i * s_i) * G
#       = product(k_i * G) * (sum(L_i * s_i) * G)^c
#       = product(R_i) * (secret_key * G)^c  (since sum(L_i * s_i) reconstructs the secret)
#       = R * (Y^c)
