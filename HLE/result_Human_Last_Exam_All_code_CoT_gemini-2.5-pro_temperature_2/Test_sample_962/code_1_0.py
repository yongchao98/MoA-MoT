import random
import math

# Use a large prime number to simulate the finite field for our cryptographic operations
PRIME = 428991448375630485085189344465493393247

def poly_eval(poly, x):
    """Evaluates a polynomial f(x) at a given point x."""
    res = 0
    for i, coeff in enumerate(poly):
        res = (res + coeff * pow(x, i, PRIME)) % PRIME
    return res

def extended_gcd(a, b):
    """Extended Euclidean Algorithm to find modular inverse."""
    if a == 0:
        return b, 0, 1
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y

def mod_inverse(n, prime):
    """Calculates the modular inverse of n modulo prime."""
    _, inv, _ = extended_gcd(n, prime)
    return inv % prime

def reconstruct_secret(shares, participating_ids):
    """
    Reconstructs the secret from a list of shares using Lagrange Interpolation.
    The secret corresponds to the polynomial's value at x=0.
    """
    if len(participating_ids) < t:
        raise ValueError("Not enough shares to reconstruct secret")

    k = len(participating_ids)
    sum_val = 0
    equation_parts = []

    print("--- Reconstructing Final Signature ---")
    print("Final Signature = Σ (Lagrange_i * Partial_Signature_i)\n")
    
    for i_idx, i in enumerate(participating_ids):
        # Calculate Lagrange coefficient for participant i
        numerator = 1
        denominator = 1
        for j_idx, j in enumerate(participating_ids):
            if i_idx != j_idx:
                numerator = (numerator * (0 - j)) % PRIME
                denominator = (denominator * (i - j)) % PRIME
        
        # L_i = numerator / denominator
        lagrange_coeff = (numerator * mod_inverse(denominator, PRIME)) % PRIME
        partial_sig = shares[i]
        
        # Add to the sum
        term = (lagrange_coeff * partial_sig) % PRIME
        sum_val = (sum_val + term) % PRIME

        # Print the equation part for this participant
        print(f"Term for Participant {i}:")
        print(f"  Partial Signature: {partial_sig}")
        print(f"  Lagrange Coefficient L_{i}: {lagrange_coeff}")
        print(f"  Term Value: ({lagrange_coeff} * {partial_sig}) mod PRIME = {term}")
        equation_parts.append(f"({lagrange_coeff} * {partial_sig})")

    # Display the final combined equation
    final_equation_str = " + ".join(equation_parts)
    print("\nFinal Equation:")
    print(f"Reconstructed Signature = ( {final_equation_str} ) mod {PRIME}")

    return sum_val

# --- Scheme Parameters ---
# n: Total number of participants
# t: Threshold of participants required to sign
n = 10
t = 4
print(f"System Parameters: n={n}, t={t}\n")

# --- 1. Distributed Key Generation (DKG) Simulation ---
# A dealer (or a DKG protocol) creates a secret polynomial of degree t-1.
# The master secret key `s` is the constant term f(0).
secret_polynomial = [random.randrange(1, PRIME) for _ in range(t)]
master_secret_key = secret_polynomial[0]

print(f"--- Key Generation ---")
print(f"Generated a secret polynomial of degree {t-1}.")
print(f"Master Secret Key (s = f(0)): {master_secret_key}\n")

# Generate secret shares for all n participants.
# Participant i gets share s_i = f(i). Participant IDs are 1, 2, 3, ...
secret_shares = {i: poly_eval(secret_polynomial, i) for i in range(1, n + 1)}
print("Generated secret shares for all participants.\n")

# --- 2. Signing Simulation ---
# Let's say a group of `t` participants decide to sign a message.
participating_ids = random.sample(range(1, n + 1), t)
participating_ids.sort()
print(f"--- Signing Process ---")
print(f"A group of {t} participants will sign: {participating_ids}\n")

# A message hash is simulated as a random integer. In a real system,
# this would be the output of a cryptographic hash function like SHA-256
# mapped to a point on an elliptic curve, then represented as a number.
message_hash = random.randrange(1, PRIME)
print(f"Message Hash (H(m)): {message_hash}\n")

# Each participant computes their partial signature.
# In a real BLS scheme, this is σ_i = s_i * H(m), where H(m) is a point and
# * is scalar multiplication. We simulate this with modular multiplication.
partial_signatures = {}
for pid in participating_ids:
    s_i = secret_shares[pid]
    partial_signatures[pid] = (s_i * message_hash) % PRIME

# --- 3. Signature Aggregation ---
# The aggregator collects the partial signatures and reconstructs the full signature.
# Note: For simplicity, we are passing the partial signatures indexed by participant ID.
reconstructed_signature = reconstruct_secret(partial_signatures, participating_ids)

# --- 4. Verification ---
# Anyone can compute the expected signature using the master public key.
# Expected Signature = master_secret_key * message_hash
# In a real system: verify(PK, σ, m) -> e(σ, G) == e(H(m), PK)
expected_signature = (master_secret_key * message_hash) % PRIME

print("\n--- Verification ---")
print(f"Expected Final Signature (s * H(m)): {expected_signature}")
print(f"Reconstructed Final Signature:          {reconstructed_signature}")

if reconstructed_signature == expected_signature:
    print("\nSUCCESS: The reconstructed signature is valid!")
else:
    print("\nFAILURE: The reconstructed signature is invalid!")