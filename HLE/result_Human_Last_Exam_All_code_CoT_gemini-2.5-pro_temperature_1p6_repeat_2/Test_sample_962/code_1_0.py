import random

# WARNING: This code is for educational purposes only.
# It demonstrates the mathematical concepts of a threshold signature scheme
# but is NOT cryptographically secure. Do not use in production.

# --- Parameters ---
# Prime modulus for our finite field Z_p
p = 65537 # A small Fermat prime for demonstration
# Generator g of the multiplicative group of Z_p
g = 3
# Total number of participants
n = 5
# Threshold of participants required to sign
t = 3

print(f"Scheme Parameters:\n  p (prime modulus) = {p}\n  g (generator) = {g}\n  n (total participants) = {n}\n  t (threshold) = {t}\n")

# --- Helper Functions ---
def power(base, exp, mod):
    """Computes (base^exp) % mod efficiently."""
    return pow(base, exp, mod)

def inv(a, mod):
    """Computes the modular multiplicative inverse of a modulo mod."""
    return pow(a, -1, mod)

def lagrange_coeff(S, i, field_order):
    """
    Computes the Lagrange coefficient L_i for participant i in signer set S.
    L_i = product_{j in S, j != i} (j / (j - i))
    Note: Participants are identified by indices 1, 2, 3...
    The field order is p-1 for exponents.
    """
    num = 1
    den = 1
    for j in S:
        if i == j:
            continue
        num = (num * j) % field_order
        den = (den * (j - i)) % field_order
    return (num * inv(den, field_order)) % field_order

# A simple, insecure "hash" function for demonstration.
# In a real system, this would be a secure cryptographic hash like SHA-256.
def pseudo_hash(*args):
    """Combine multiple numbers into one. Not a real hash function!"""
    s = sum(int(x) for x in args)
    return s % (p - 1)

# --- 1. Key Generation (Performed by a trusted dealer) ---
print("--- 1. Key Generation Phase ---")
# The secret polynomial f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
# The main secret is a_0 = f(0)
coefficients = [random.randint(1, p - 2) for _ in range(t - 1)]
secret_key = random.randint(1, p-2)
coefficients.insert(0, secret_key) # f(0) is the secret key

print(f"Secret polynomial coefficients: {coefficients}")
print(f"Master Secret Key (f(0)): {secret_key}")

def f(x):
    """Evaluates the secret polynomial at x."""
    y = 0
    # Evaluate using Horner's method for efficiency
    for coeff in reversed(coefficients):
        y = (y * x + coeff) % (p - 1)
    return y

# Generate private key shares and public verification keys
private_shares = {}
public_verification_keys = {}
for i in range(1, n + 1):
    private_shares[i] = f(i)
    public_verification_keys[i] = power(g, private_shares[i], p)

# The group's public key
public_key = power(g, secret_key, p)
print(f"Group Public Key (Y = g^secret): {public_key}\n")

# --- 2. Signing Protocol ---
print("--- 2. Signing Phase ---")

# The set of participants who will sign the message
signers = list(range(1, t + 1))
print(f"Participants signing: {signers}")

# The message to be signed
message = 12345
print(f"Message to sign (m): {message}\n")


# === ROUND 1: Commitment ===
print("--- ROUND 1: Commitments ---")
nonces = {} # To store (d_i, e_i)
commitments = {} # To store (D_i, E_i)

for i in signers:
    # Each signer generates two private random nonces
    d_i = random.randint(1, p - 2)
    e_i = random.randint(1, p - 2)
    nonces[i] = (d_i, e_i)

    # And computes the corresponding public commitments
    D_i = power(g, d_i, p)
    E_i = power(g, e_i, p)
    commitments[i] = (D_i, E_i)
    print(f"  P{i}: Generated nonces (d_{i}, e_{i}) and commitments (D_{i}, E_{i}) = ({D_i}, {E_i})")
print("")

# === ROUND 2: Signing ===
print("--- ROUND 2: Aggregation and Response ---")

# The signing coordinator gathers all commitments
commitment_list = []
for i in signers:
    commitment_list.extend(commitments[i])

# Step 1: Coordinator computes binding values for each participant
binding_values = {}
for i in signers:
    # The binding value rho_i = H(i, m, {D_j, E_j})
    # This prevents malleability attacks.
    binding_values[i] = pseudo_hash(i, message, *commitment_list)
print(f"Coordinator computed binding values (rho): {binding_values}")

# Step 2: Coordinator computes the group commitment R
# R = Product_{i in S} (D_i * E_i^{rho_i})
R = 1
for i in signers:
    D_i, E_i = commitments[i]
    rho_i = binding_values[i]
    term = (D_i * power(E_i, rho_i, p)) % p
    R = (R * term) % p
print(f"Coordinator computed group commitment (R): {R}")

# Step 3: Coordinator computes the challenge c = H(Y, R, m)
challenge = pseudo_hash(public_key, R, message)
print(f"Coordinator computed challenge (c): {challenge}\n")


# Step 4: Each participant computes their partial signature z_i
# z_i = d_i + e_i * rho_i + c * L_i * s_i mod (p-1)
partial_signatures = {}
for i in signers:
    d_i, e_i = nonces[i]
    rho_i = binding_values[i]
    s_i = private_shares[i]

    # Calculate this participant's Lagrange coefficient
    lambda_i = lagrange_coeff(signers, i, p - 1)

    # Compute the response
    z_i = (d_i + (e_i * rho_i) + (challenge * lambda_i * s_i)) % (p - 1)
    partial_signatures[i] = z_i
    print(f"  P{i}: Computed Lagrange coeff L_{i}={lambda_i}, partial signature z_{i}={z_i}")

# Step 5: Coordinator aggregates the partial signatures
# z = sum(z_i)
z = sum(partial_signatures.values()) % (p-1)
print("\nCoordinator aggregates partial signatures.")

print("\n--- 3. Final Signature ---")
print(f"The final aggregated signature is (R, z):")
print(f"  R = {R}")
print(f"  z = {z}")

# --- 4. Verification ---
print("\n--- 4. Verification Step ---")
print("A verifier checks if g^z == R * Y^c mod p")

# Left-hand side of the equation: g^z
lhs = power(g, z, p)

# Right-hand side of the equation: R * Y^c
rhs = (R * power(public_key, challenge, p)) % p

print(f"\nFinal Equation for Verification:")
print(f"Is {g}^{z} mod {p}  ==  ({R} * {public_key}^{challenge}) mod {p} ?")
print(f"Calculating LHS (g^z): {lhs}")
print(f"Calculating RHS (R * Y^c): {rhs}")

if lhs == rhs:
    print("\nResult: SUCCESS! The signature is valid.")
    <<<SUCCESS>>>
else:
    print("\nResult: FAILURE! The signature is invalid.")
    <<<FAILURE>>>