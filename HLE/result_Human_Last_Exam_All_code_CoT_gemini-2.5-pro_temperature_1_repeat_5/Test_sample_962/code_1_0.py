import random
import hashlib

# A simple class to simulate a finite field for demonstration
class Field:
    def __init__(self, p):
        self.p = p
    def add(self, a, b):
        return (a + b) % self.p
    def mul(self, a, b):
        return (a * b) % self.p
    def inv(self, a):
        return pow(a, -1, self.p)
    def sub(self, a, b):
        return (a - b + self.p) % self.p

# --- Protocol Parameters (Illustrative) ---
# Prime modulus for a finite field (a subgroup order in a real implementation)
q = 251
F = Field(q)
# A generator element
g = 3
# Threshold parameters
t = 3
n = 5

# --- Helper Functions ---
def H(data):
    """A simple hash function for the demonstration."""
    # In a real implementation, this would be a cryptographic hash like SHA-256
    # and its output would be mapped to an integer mod q.
    if not isinstance(data, bytes):
        data = str(data).encode()
    return int(hashlib.sha256(data).hexdigest(), 16) % q

def get_lagrange_coeff(S, j):
    """Calculates the Lagrange coefficient L_j for participant j in set S."""
    S_list = list(S)
    num, den = 1, 1
    for i in S_list:
        if i != j:
            num = F.mul(num, i)
            den = F.mul(den, F.sub(i, j))
    return F.mul(num, F.inv(den))

# --- 1. Key Generation (Simulated by a Trusted Dealer) ---
print("--- 1. Key Generation ---")
# The master secret key
master_secret_key = random.randint(1, q - 1)
print(f"Master Secret Key (x): {master_secret_key}")

# A (t-1) degree polynomial f(z) = a_0 + a_1*z + ... + a_{t-1}*z^{t-1}
# where a_0 is the master secret key
coeffs = [master_secret_key] + [random.randint(1, q - 1) for _ in range(t - 1)]

def poly(z):
    res = 0
    for i in range(t):
        res = F.add(res, F.mul(coeffs[i], pow(z, i, q)))
    return res

# Generate secret shares for n participants
secret_shares = {i: poly(i) for i in range(1, n + 1)}
public_verif_keys = {i: pow(g, s, q) for i, s in secret_shares.items()}
group_public_key = pow(g, master_secret_key, q)

print(f"Group Public Key (Y = g^x): {group_public_key}")
print(f"Secret Shares (s_i): {secret_shares}")
print("-" * 20 + "\n")

# --- 2. Signing Protocol ---
print("--- 2. Signing Protocol ---")
message = "This is a test message"
print(f"Message (m): '{message}'")

# Assume participants {1, 2, 3} are signing (a set of size t)
S = {1, 2, 3}
print(f"Signing Participants (S): {S}")

# --- Round 1: Commitment ---
print("\n--- Round 1: Commitment ---")
nonces_d = {i: random.randint(1, q - 1) for i in S}
nonces_e = {i: random.randint(1, q - 1) for i in S}

commits_D = {i: pow(g, d, q) for i, d in nonces_d.items()}
commits_E = {i: pow(g, e, q) for i, e in nonces_e.items()}

print("Each participant i in S generates random nonces (d_i, e_i) and broadcasts commitments (D_i, E_i).")
for i in S:
    print(f"  P_{i}: d_{i}={nonces_d[i]}, e_{i}={nonces_e[i]} -> D_{i}={commits_D[i]}, E_{i}={commits_E[i]}")

# --- Round 2: Signing ---
print("\n--- Round 2: Signing ---")

# Step 2a: Compute binding values (rho_i) and group commitment (R)
print("\nStep 2a: Compute binding values and group commitment")
# In FROST, rho is derived from the message and all commitments
all_commits_list = [v for k in S for v in (commits_D[k], commits_E[k])]
binding_base_hash = H(str((list(S), message, all_commits_list)))
rho = {i: H(str((binding_base_hash, i))) for i in S}

# Calculate group commitment R
R = 1
r_equation_parts = []
for i in S:
    E_rho = pow(commits_E[i], rho[i], q)
    R = F.mul(R, F.mul(commits_D[i], E_rho))
    r_equation_parts.append(f"(D_{i} * E_{i}^rho_{i})")

print(f"Binding values (rho_i): {rho}")
print(f"Group Commitment R = {' * '.join(r_equation_parts)} (mod q)")
print(f"R = {R}")

# Step 2b: Compute challenge (c)
c = H(str((R, group_public_key, message)))
print(f"\nStep 2b: Compute Challenge c = H(R, Y, m) = {c}")

# Step 2c: Each participant computes their response (z_i)
print("\nStep 2c: Each participant computes their signature share z_i")
print("z_i = d_i + (e_i * rho_i) + c * L_i * s_i (mod q)")

lagrange_coeffs = {i: get_lagrange_coeff(S, i) for i in S}
z_shares = {}
for i in S:
    s_i = secret_shares[i]
    L_i = lagrange_coeffs[i]
    term1 = nonces_d[i]
    term2 = F.mul(nonces_e[i], rho[i])
    term3 = F.mul(F.mul(c, L_i), s_i)
    z_i = F.add(term1, F.add(term2, term3))
    z_shares[i] = z_i
    print(f"  P_{i}: L_{i}={L_i}, z_{i} = {term1} + {term2} + {term3} = {z_i}")

# Step 2d: Aggregator combines shares into the final signature
print("\nStep 2d: Aggregator assembles the final signature (R, s)")
s = sum(z_shares.values()) % q
final_signature = (R, s)

z_share_names = [f"z_{i}" for i in S]
print(f"Final signature s = {' + '.join(z_share_names)} (mod q)")
print(f"s = {' + '.join(map(str, z_shares.values()))} mod {q} = {s}")
print(f"\nFinal Signature (R, s) = {final_signature}")
print("-" * 20 + "\n")

# --- 3. Verification ---
print("--- 3. Verification ---")
print("Verifier checks if g^s == R * Y^c (mod q)")
R_val, s_val = final_signature
lhs = pow(g, s_val, q)
rhs = F.mul(R_val, pow(group_public_key, c, q))

print(f"g^s = {g}^{s_val} mod {q} = {lhs}")
print(f"R * Y^c = {R_val} * {group_public_key}^{c} mod {q} = {rhs}")

if lhs == rhs:
    print("\nSUCCESS: The signature is valid.")
else:
    print("\nFAILURE: The signature is invalid.")
