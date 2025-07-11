import random
import hashlib

# Step 1: Setup
# We use a Schnorr-like scheme over integers modulo a prime p.
# q is a prime that divides p-1.
# g is a generator of a subgroup of order q.
# For this example, p=23, q=11, and g=4. (4^11 mod 23 = 1)
p = 23
q = 11
g = 4
n = 5  # Total number of participants
t = 3  # Threshold of participants required to sign

def H(data):
    """A simple hash function for the simulation."""
    sha = hashlib.sha256()
    if isinstance(data, (list, tuple)):
        for item in sorted(data):
            sha.update(str(item).encode())
    else:
        sha.update(str(data).encode())
    return int(sha.hexdigest(), 16) % q

def power(base, exp, mod):
    """Modular exponentiation: (base^exp) % mod"""
    return pow(base, exp, mod)

def inverse(a, mod):
    """Modular inverse: (a^-1) % mod"""
    return pow(a, -1, mod)

def generate_polynomial(degree, intercept):
    """Generates a random polynomial f(x) of a given degree where f(0) = intercept."""
    coefficients = [intercept] + [random.randint(1, q - 1) for _ in range(degree)]
    return coefficients

def evaluate_polynomial(coeffs, x):
    """Evaluates a polynomial at a given point x."""
    res = 0
    for i, coeff in enumerate(coeffs):
        res = (res + coeff * power(x, i, q)) % q
    return res

def lagrange_coefficient(i, signers):
    """Calculates the Lagrange basis polynomial coefficient lambda_i for participant i."""
    num, den = 1, 1
    for j in signers:
        if i != j:
            num = (num * j) % q
            den = (den * (j - i)) % q
    return (num * inverse(den, q)) % q

class Participant:
    def __init__(self, id, secret_share):
        self.id = id
        self.sk_i = secret_share
        self.d_i = None  # Nonce 1
        self.e_i = None  # Nonce 2
        self.D_i = None  # Commitment to d_i
        self.E_i = None  # Commitment to e_i

    def __repr__(self):
        return f"Participant(id={self.id})"

    def generate_nonces_and_commitments(self):
        """Round 1: Generate nonces and their commitments."""
        self.d_i = random.randint(1, q - 1)
        self.e_i = random.randint(1, q - 1)
        self.D_i = power(g, self.d_i, p)
        self.E_i = power(g, self.e_i, p)
        return self.D_i, self.E_i

    def generate_signature_share(self, message, commitments, group_pk, signers):
        """Round 2: Generate the partial signature z_i."""
        # Calculate binding value b
        commitment_list = [c for c_id in sorted(commitments.keys()) for c in commitments[c_id]]
        b = H((group_pk, message, commitment_list))

        # Calculate group commitment R
        R_val = 1
        for j in signers:
            D_j, E_j = commitments[j]
            R_val = (R_val * D_j * power(E_j, b, p)) % p

        # Calculate challenge c
        c = H((group_pk, R_val, message))

        # Calculate partial signature z_i
        lambda_i = lagrange_coefficient(self.id, signers)
        s_i_c_lambda = (c * self.sk_i * lambda_i) % q
        r_i = (self.d_i + self.e_i * b) % q
        z_i = (r_i + s_i_c_lambda) % q
        
        print(f"  Participant {self.id}:")
        print(f"    - Nonces (d, e): ({self.d_i}, {self.e_i})")
        print(f"    - Binding value b: {b}")
        print(f"    - Lagrange coeff lambda_{self.id}: {lambda_i}")
        print(f"    - Partial signature z_{self.id} = (d_i + e_i*b) + c*sk_i*lambda_i = {z_i}")
        
        return z_i

# --- Main Simulation ---

print("--- 1. Key Generation ---")
# A trusted dealer generates the master secret and splits it.
master_sk = random.randint(1, q - 1)
master_pk = power(g, master_sk, p)

print(f"Global Parameters: p={p}, q={q}, g={g}")
print(f"Threshold Scheme: t={t}, n={n}")
print(f"Master Secret Key (f(0)): {master_sk}")
print(f"Master Public Key (g^f(0)): {master_pk}\n")

poly_coeffs = generate_polynomial(t - 1, master_sk)
participants = []
for i in range(1, n + 1):
    secret_share = evaluate_polynomial(poly_coeffs, i)
    participants.append(Participant(i, secret_share))
    print(f"Participant {i} gets secret share f({i}): {secret_share}")

print("\n--- 2. Signing Protocol ---")
message_to_sign = "hello world"
# A random subset of t participants will sign
signing_participants = random.sample(participants, t)
signer_ids = sorted([p.id for p in signing_participants])
print(f"Message: '{message_to_sign}'")
print(f"A group of {t} participants will sign: IDs {signer_ids}\n")

# --- Round 1: Commitment ---
print("--- Round 1: Commitment Generation (Preprocessing) ---")
commitments = {}
for p_signer in signing_participants:
    D_i, E_i = p_signer.generate_nonces_and_commitments()
    commitments[p_signer.id] = (D_i, E_i)
    print(f"Participant {p_signer.id} broadcasts commitments (D_{p_signer.id}, E_{p_signer.id}): ({D_i}, {E_i})")

# --- Round 2: Signature Share Generation ---
print("\n--- Round 2: Signature Share Generation ---")
partial_signatures = {}

for p_signer in signing_participants:
    z_i = p_signer.generate_signature_share(message_to_sign, commitments, master_pk, signer_ids)
    partial_signatures[p_signer.id] = z_i

# --- 3. Aggregation and Verification ---
print("\n--- 3. Aggregation and Verification ---")

# The aggregator (can be one of the signers) collects partial signatures
s = sum(partial_signatures.values()) % q
print(f"Aggregator computes final signature 's' = sum(z_i) mod q: {s}")

# The aggregator also computes the group commitment R
commitment_list = [c for c_id in sorted(commitments.keys()) for c in commitments[c_id]]
binding_val = H((master_pk, message_to_sign, commitment_list))
R = 1
for signer_id in signer_ids:
    D_j, E_j = commitments[signer_id]
    R = (R * D_j * power(E_j, binding_val, p)) % p
print(f"Aggregator computes group commitment 'R': {R}")

# The challenge c is re-calculated for verification
challenge = H((master_pk, R, message_to_sign))
print(f"Aggregator computes challenge 'c' = H(PK, R, m): {challenge}")

# The final signature is the tuple (R, s)
print(f"\nFinal Signature is (R, s) = ({R}, {s})")

# Verification
print("\n--- Final Verification Equation: g^s == R * PK^c (mod p) ---")
lhs = power(g, s, p)
rhs_term2 = power(master_pk, challenge, p)
rhs = (R * rhs_term2) % p

print(f"Calculating left side:  g^s = {g}^{s} = {lhs} (mod {p})")
print(f"Calculating right side: R * PK^c = {R} * {master_pk}^{challenge} = {rhs} (mod {p})")

is_valid = lhs == rhs
print(f"\nSignature is valid: {is_valid}")

final_tuple = (R,s)
<<<final_tuple>>>