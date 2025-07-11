import random
import hashlib
from math import gcd

# --- Helper Functions for Modular Arithmetic ---
# (These are needed for Shamir's Secret Sharing and Schnorr signatures)

def extended_gcd(a, b):
    """Returns (gcd, x, y) such that a*x + b*y = gcd."""
    if a == 0:
        return (b, 0, 1)
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return (d, x, y)

def mod_inverse(k, p):
    """Returns the modular inverse of k mod p."""
    if gcd(k, p) != 1:
        return None
    _, x, _ = extended_gcd(k, p)
    return x % p

def lagrange_interpolate(x, points, p):
    """
    Finds the y-value for a given x using Lagrange interpolation over a prime field.
    'points' is a list of (x, y) coordinates.
    """
    x_coords, y_coords = zip(*points)
    L = 0
    for i in range(len(points)):
        num, den = 1, 1
        for j in range(len(points)):
            if i == j:
                continue
            num = (num * (x - x_coords[j])) % p
            den = (den * (x_coords[i] - x_coords[j])) % p
        term = (y_coords[i] * num * mod_inverse(den, p)) % p
        L = (L + term) % p
    return L

# --- Simulation of a 2-Round Threshold Signature Scheme (FROST-inspired) ---

def run_threshold_signature_demo():
    """
    Demonstrates a simplified 2-round t-out-of-n threshold signature scheme.
    """
    print("--- 1. Setup Phase ---")
    # n: Total number of parties in the key-holding group
    n = 5
    # t: Threshold of signers required to produce a signature
    t = 3
    # A large prime for the finite field (in real crypto, this would be much larger)
    p = 65537
    # A generator g for the group
    g = 3

    print(f"Total parties (n): {n}, Threshold (t): {t}\n")

    # --- 2. Key Generation Phase (Simulated with a Trusted Dealer) ---
    print("--- 2. Key Generation Phase (Trusted Dealer) ---")
    # In a real system, a Distributed Key Generation (DKG) protocol is used.
    # The dealer creates a secret polynomial of degree t-1.
    # f(x) = a_0 + a_1*x + ... + a_{t-1}*x^{t-1}
    # The master secret key is sk = a_0 = f(0)
    coeffs = [random.randint(1, p - 1) for _ in range(t)]
    secret_key = coeffs[0]

    def polynomial(x):
        val = 0
        for i, c in enumerate(coeffs):
            val = (val + c * (x**i)) % p
        return val

    # Generate secret shares for each of the n parties
    secret_shares = {i: polynomial(i) for i in range(1, n + 1)}
    print(f"Secret shares distributed to {n} parties.")
    # for i, s in secret_shares.items():
    #     print(f"  - Party {i} receives secret share s_{i} = {s}")

    # The single group public key is Y = g^secret_key
    public_key = pow(g, secret_key, p)
    print(f"Group Public Key (Y): {public_key}\n")

    # --- 3. Signing Phase ---
    print("--- 3. Signing Phase ---")
    # A message to be signed
    message = "This is a test message for threshold signing"
    
    # A subset of t parties will participate in signing
    signer_indices = random.sample(range(1, n + 1), t)
    print(f"A group of {t} parties will sign: {signer_indices}")

    # -- Round 1: Commitments --
    print("\n-- Round 1: Commitments --")
    nonces = {}
    commitments = {}
    for i in signer_indices:
        # Each signer generates a secret nonce k_i
        k_i = random.randint(1, p - 1)
        nonces[i] = k_i
        # Each signer computes and "broadcasts" a public commitment R_i = g^k_i
        R_i = pow(g, k_i, p)
        commitments[i] = R_i
        print(f"  - Party {i} generates nonce k_{i} and broadcasts commitment R_{i} = {R_i}")

    # -- Round 2: Signature Shares --
    print("\n-- Round 2: Signature Shares --")
    # All signers have now received the commitments from Round 1.
    
    # Each signer computes the group commitment R
    group_commitment_R = 1
    for R_i in commitments.values():
        group_commitment_R = (group_commitment_R * R_i) % p
    print(f"All signers compute the same Group Commitment (R): {group_commitment_R}")

    # Each signer computes the challenge c = H(R, Y, m)
    h = hashlib.sha256()
    h.update(str(group_commitment_R).encode())
    h.update(str(public_key).encode())
    h.update(message.encode())
    challenge_c = int(h.hexdigest(), 16) % p
    print(f"All signers compute the same Challenge (c): {challenge_c}")

    partial_signatures = {}
    for i in signer_indices:
        # Each signer computes their Lagrange coefficient for this group of signers
        # This allows reconstruction of the master secret at x=0 from the shares
        points = [(j, 0) for j in signer_indices] # Use dummy y-values
        l_i = lagrange_interpolate(0, [(j, 1 if j==i else 0) for j in signer_indices], p)

        # Each signer computes their partial signature: z_i = k_i + c * l_i * s_i
        s_i = secret_shares[i]
        k_i = nonces[i]
        z_i = (k_i + challenge_c * l_i * s_i) % p
        partial_signatures[i] = z_i
        print(f"  - Party {i} computes and broadcasts partial signature z_{i} = {z_i}")

    # --- 4. Signature Aggregation ---
    print("\n--- 4. Signature Aggregation ---")
    # An aggregator collects the partial signatures and sums them up
    final_signature_z = sum(partial_signatures.values()) % p
    
    print("Aggregator sums partial signatures to get the final signature component z.")
    print(f"Final Signature: (R = {group_commitment_R}, z = {final_signature_z})")

    # --- 5. Verification ---
    print("\n--- 5. Verification ---")
    print("A verifier checks if g^z == R * (Y^c) mod p")
    
    # Left side of the equation: g^z
    left_side = pow(g, final_signature_z, p)
    # Right side of the equation: R * (Y^c)
    right_side = (group_commitment_R * pow(public_key, challenge_c, p)) % p

    print("\nFinal Equation Check:")
    print(f"  g^z mod p:   {g}^{final_signature_z} mod {p} = {left_side}")
    print(f"  R*Y^c mod p: {group_commitment_R} * {public_key}^{challenge_c} mod {p} = {right_side}")
    
    if left_side == right_side:
        print("\nResult: SUCCESS! The signature is valid.")
    else:
        print("\nResult: FAILURE! The signature is invalid.")

# Run the demonstration
run_threshold_signature_demo()
