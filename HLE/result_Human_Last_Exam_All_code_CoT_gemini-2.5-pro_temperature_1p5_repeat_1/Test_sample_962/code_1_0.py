import hashlib
import random

# --- Cryptographic Primitives (using modular arithmetic for simplicity) ---
# In a real implementation, these would be elliptic curve operations.

# Prime modulus and generator for our finite field
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
G = 2

def sha256_as_int(data_str):
    """Hashes a string and returns it as an integer."""
    return int(hashlib.sha256(data_str.encode()).hexdigest(), 16)

# Simulates scalar multiplication: k * BasePoint -> g^k mod p
def scalar_mult(k):
    return pow(G, k, P)

# Simulates point addition: P1 + P2 -> v1 * v2 mod p
def point_add(v1, v2):
    return (v1 * v2) % P

# --- Shamir's Secret Sharing Polynomial ---
class Polynomial:
    def __init__(self, degree, secret=None):
        self.degree = degree
        # Generate random coefficients. The secret is the constant term f(0).
        self.coeffs = [random.randint(1, P - 1) for _ in range(degree)]
        if secret:
            self.coeffs.insert(0, secret)
        else:
            self.coeffs.insert(0, random.randint(1, P - 1))

    def evaluate(self, x):
        """Evaluate the polynomial at a given point x."""
        res = 0
        for i in range(self.degree + 1):
            res = (res + self.coeffs[i] * pow(x, i, P)) % P
        return res

# --- Lagrange Interpolation Coefficient ---
def lagrange_coeff(i, S):
    """Computes the Lagrange coefficient lambda_i for participant i in signing set S."""
    num = 1
    den = 1
    for j in S:
        if i == j:
            continue
        num = (num * j) % P
        den = (den * (j - i)) % P
    return (num * pow(den, -1, P)) % P


# --- Simulation Setup ---
N = 5  # Total number of participants
T = 3  # Threshold required to sign

print(f"--- System Setup: {T}-out-of-{N} Threshold Signature ---")

# 1. Trusted Dealer generates keys (replaces a DKG)
master_secret_key = random.randint(1, P - 1)
poly = Polynomial(T - 1, master_secret_key)
master_public_key = scalar_mult(master_secret_key)

print(f"Master Public Key (Y): {master_public_key}\n")

secret_shares = {i: poly.evaluate(i) for i in range(1, N + 1)}
public_shares = {i: scalar_mult(s) for i, s in secret_shares.items()}

# 2. A signing set of size T is chosen
signing_participants_ids = random.sample(range(1, N + 1), T)
print(f"--- Signing Initiated by Participants {signing_participants_ids} ---")
message = "This is a very important message to sign"
print(f"Message (m): '{message}'\n")

# --- Signing Round 1: Commitments ---
print("--- Round 1: Nonce Commitments ---")
# Each participant generates secret nonces (d, e) and public commitments (D, E)
nonces = {} # Stores (d_i, e_i) for each participant i
commitments = {} # Stores (D_i, E_i) for each participant i

for i in signing_participants_ids:
    d_i = random.randint(1, P - 1)
    e_i = random.randint(1, P - 1)
    nonces[i] = (d_i, e_i)
    
    D_i = scalar_mult(d_i)
    E_i = scalar_mult(e_i)
    commitments[i] = (D_i, E_i)
    print(f"Participant {i} broadcasts commitments (D_{i}, E_{i}): ({D_i}, {E_i})")

# --- Signing Round 2: Partial Signatures ---
print("\n--- Round 2: Partial Signature Generation ---")

# Each participant computes the group commitment R
# First, create the binding value rho_i for each participant
# rho_i = H(i, m, {D_j, E_j}_j)
binding_factors = {}
commitment_list_str = "".join([str(v) for c in commitments.values() for v in c])

for i in signing_participants_ids:
    hash_input = f"{i}{message}{commitment_list_str}"
    binding_factors[i] = sha256_as_int(hash_input)

# R = sum(D_j + rho_j * E_j) -> Product(D_j * E_j^rho_j)
group_commitment_R = 1
for j in signing_participants_ids:
    D_j, E_j = commitments[j]
    rho_j = binding_factors[j]
    term = (D_j * pow(E_j, rho_j, P)) % P
    group_commitment_R = point_add(group_commitment_R, term)

print(f"Group Commitment (R): {group_commitment_R}")

# Each participant computes the challenge c = H(Y, R, m)
challenge_hash_input = f"{master_public_key}{group_commitment_R}{message}"
challenge_c = sha256_as_int(challenge_hash_input)
print(f"Challenge (c): {challenge_c}\n")

# Each participant computes and broadcasts their partial signature z_i
# z_i = d_i + e_i * rho_i + c * lambda_i * s_i
partial_signatures = {}
for i in signing_participants_ids:
    d_i, e_i = nonces[i]
    rho_i = binding_factors[i]
    lambda_i = lagrange_coeff(i, signing_participants_ids)
    s_i = secret_shares[i]

    z_i = (d_i + e_i * rho_i + challenge_c * lambda_i * s_i) % P
    partial_signatures[i] = z_i
    print(f"Participant {i} computes partial signature z_{i}: {z_i}")


# --- Aggregation and Verification ---
print("\n--- Aggregation and Verification ---")

# An aggregator sums the partial signatures to get the final signature z
final_z = sum(partial_signatures.values()) % P
print(f"Aggregated Signature (z): {final_z}")
print(f"Final Signature is the pair (R, z): ({group_commitment_R}, {final_z})")

# Verification: Check if z * G == R + c * Y
# In our modular arithmetic simulation: g^z == R * Y^c mod p

# Left side of the equation: z * G -> g^z
lhs = scalar_mult(final_z)

# Right side of the equation: R + c*Y -> R * Y^c
rhs = point_add(group_commitment_R, pow(master_public_key, challenge_c, P))

print("\n--- Final Verification Equation ---")
print(f"Does g^z mod p == (R * Y^c) mod p ?")
print(f"Left side (g^z):  {lhs}")
print(f"Right side (R*Y^c): {rhs}")

if lhs == rhs:
    print("\nSUCCESS: The signature is valid!")
else:
    print("\nFAILURE: The signature is invalid!")
