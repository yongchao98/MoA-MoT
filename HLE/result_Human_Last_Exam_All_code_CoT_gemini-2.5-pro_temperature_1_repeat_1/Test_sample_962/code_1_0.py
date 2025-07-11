import random
import hashlib

# Step 1: Setup
# We simulate a cryptographic group using modular arithmetic.
# In a real implementation, this would be an elliptic curve group with a pairing.
# q is the order of the group, p is the prime field modulus.
q = 115792089237316195423570985008687907853269984665640564039457584007908834671663
p = 115792089237316195423570985008687907853269984665640564039457584007913129639935
g = 2 # A generator of the group

# n-out-of-t threshold scheme parameters
n = 5 # Total number of parties
t = 3 # Threshold required to sign

def extended_gcd(a, b):
    """Extended Euclidean Algorithm to find modular inverse."""
    if a == 0:
        return b, 0, 1
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y

def modInverse(a, m):
    """Modular inverse of a mod m."""
    d, x, y = extended_gcd(a, m)
    if d != 1:
        raise Exception('modular inverse does not exist')
    return x % m

def hash_to_int(message):
    """A simple 'hash-to-group' simulation."""
    return int.from_bytes(hashlib.sha256(message.encode()).digest(), 'big') % q

# Step 2: Distributed Key Generation (Simulated by a Dealer)
print(f"--- Setting up a {t}-out-of-{n} Threshold Scheme ---\n")

# The dealer creates a secret polynomial of degree t-1
# The secret key is the constant term, polynomial[0]
polynomial = [random.randint(1, q-1) for _ in range(t)]
secret_key = polynomial[0]
print(f"Dealer's secret polynomial coefficients (secret key is the first): {polynomial}")

def evaluate_polynomial(poly, x):
    """Evaluates polynomial f(x) at a given point x."""
    y = 0
    # Horner's method for polynomial evaluation
    for coeff in reversed(poly):
        y = (y * x + coeff) % q
    return y

# Generate secret shares for each of the n parties
# Each party i gets a secret share sk_i = f(i)
secret_shares = {i: evaluate_polynomial(polynomial, i) for i in range(1, n + 1)}
print("\n--- Distributing Secret Shares ---")
for i, share in secret_shares.items():
    print(f"Party {i} receives secret share: {share}")

# Calculate the public key for the group and verification keys for each share
public_key = pow(g, secret_key, p)
public_verification_keys = {i: pow(g, share, p) for i, share in secret_shares.items()}
print(f"\nGroup Public Key (g^s): {public_key}")

# Step 3: Signing Protocol
print("\n--- Initiating Signing Protocol ---")
message_to_sign = "hello world"
h = hash_to_int(message_to_sign)
print(f"Message: '{message_to_sign}'")
print(f"Hashed Message (h): {h}")

# A subset of t parties will participate in signing
# Let's choose parties 1, 3, and 5
signer_indices = [1, 3, 5]
print(f"\nRound 1: Aggregator requests signatures from parties: {signer_indices}")

# Round 2: Signers compute their partial signatures
print("\nRound 2: Signers compute and return partial signatures")
partial_signatures = {}
for i in signer_indices:
    sk_i = secret_shares[i]
    # Each signer computes sig_i = h^(sk_i)
    sig_i = pow(h, sk_i, p)
    partial_signatures[i] = sig_i
    print(f"Party {i} computes partial signature: {sig_i}")

# Step 4: Signature Aggregation
print("\n--- Aggregating Partial Signatures ---")
# The aggregator needs to compute Lagrange coefficients to combine the signatures.
# The final signature is PROD(sig_i ^ L_i) where L_i is the Lagrange coefficient.

def get_lagrange_coeff(signer_indices, i, q_mod):
    """
    Computes the Lagrange coefficient L_i(0) for a given set of signers.
    L_i(0) = product over j!=i of (j / (j-i))
    """
    numerator = 1
    denominator = 1
    for j in signer_indices:
        if i == j:
            continue
        numerator = (numerator * j) % q_mod
        denominator = (denominator * (j - i)) % q_mod
    return (numerator * modInverse(denominator, q_mod)) % q_mod

# Calculate coefficients and print them
lagrange_coeffs = {i: get_lagrange_coeff(signer_indices, i, q) for i in signer_indices}
print("Lagrange Coefficients (L_i) for reconstruction:")
for i, coeff in lagrange_coeffs.items():
    print(f"L_{i}: {coeff}")

# Combine the partial signatures
aggregated_signature = 1
print("\nFinal Equation for Aggregation: Signature = PRODUCT [ (partial_sig_i) ^ (L_i) ] mod p")
for i in signer_indices:
    partial_sig = partial_signatures[i]
    lagrange_coeff = lagrange_coeffs[i]
    # This is equivalent to sigma = PROD( (h^sk_i) ^ L_i ) = h ^ (SUM( sk_i * L_i ))
    term = pow(partial_sig, lagrange_coeff, p)
    aggregated_signature = (aggregated_signature * term) % p
    # Showing each part of the equation
    print(f"  Term for party {i}: ({partial_sig}) ^ ({lagrange_coeff}) mod p")


print(f"\nFinal Aggregated Signature: {aggregated_signature}")

# Step 5: Verification
print("\n--- Verifying the Signature ---")

# The verifier checks if the signature is valid against the group public key.
# In a real BLS scheme, this check is: e(g, signature) == e(public_key, h)
# Here we simulate it by checking if signature == h^(secret_key)
expected_signature = pow(h, secret_key, p)
print(f"Expected signature (computed with original secret): {expected_signature}")

if aggregated_signature == expected_signature:
    print("\nSUCCESS: The aggregated signature is valid!")
else:
    print("\nFAILURE: The aggregated signature is invalid.")