import random
import hashlib

# Helper function for modular inverse (to compute 1/a mod m)
def mod_inverse(a, m):
    """Computes the modular multiplicative inverse of a modulo m."""
    a = a % m
    for x in range(1, m):
        if (a * x) % m == 1:
            return x
    return 1 # Should not happen if m is prime

def get_lagrange_coeff(i, participants, q):
    """
    Calculates the Lagrange coefficient for participant i in the set of participants.
    lambda_i = product(j / (j - i)) for j in participants, j != i
    """
    numerator = 1
    denominator = 1
    for j in participants:
        if i == j:
            continue
        numerator = (numerator * j) % q
        denominator = (denominator * (j - i)) % q
    
    return (numerator * mod_inverse(denominator, q)) % q

def main():
    """Main function to simulate the t-out-of-n threshold signature scheme."""
    # --- Phase 1: Key Generation (Simulated with a Trusted Dealer) ---
    print("--- Phase 1: Key Generation ---")
    
    # 1. Public Parameters (using small numbers for demonstration)
    # In a real system, these would be very large primes.
    # p is a prime, q is a prime order of the subgroup, q divides p-1.
    p = 487 # A prime number
    q = 163 # A prime number, where 163 divides 486 (p-1)
    g = 10  # A generator of the subgroup of order q
    
    # 2. Scheme Parameters
    n = 5  # Total number of parties
    t = 3  # Threshold of parties required to sign
    
    print(f"Public Parameters: p={p}, q={q}, g={g}")
    print(f"Scheme Setup: n={n} (total parties), t={t} (threshold)\n")
    
    # 3. Trusted dealer creates a secret polynomial of degree t-1
    # f(z) = a_0 + a_1*z + ... + a_{t-1}*z^{t-1}
    # The master secret key is x = f(0) = a_0
    coefficients = [random.randint(1, q - 1) for _ in range(t)]
    master_secret_key_x = coefficients[0]
    
    def f(z):
        """The secret polynomial."""
        y = 0
        for i in range(t - 1, -1, -1):
            y = (y * z + coefficients[i]) % q
        return y

    # 4. Master public key
    master_public_key_Y = pow(g, master_secret_key_x, p)
    print(f"Master Public Key (Y): {master_public_key_Y}")
    
    # 5. Distribute shares to n parties
    # Party indices are 1, 2, ..., n
    private_shares = {i: f(i) for i in range(1, n + 1)}
    print("Private shares distributed to parties (only for demonstration):")
    for party_id, share in private_shares.items():
        print(f"  Party {party_id}: x_{party_id} = {share}")
    print("-" * 35 + "\n")

    # --- Phase 2: Signing Protocol ---
    print("--- Phase 2: Signing Protocol ---")
    message = "Hello, Threshold World!"
    print(f"Message to sign: '{message}'")
    
    # A subset of t participants will sign
    participants = random.sample(list(private_shares.keys()), t)
    participants.sort()
    print(f"Signing Participants (t={t}): {participants}\n")
    
    # -- Round 1: Commitment --
    print("-- Round 1: Commitment --")
    nonces_k = {}
    commitments_R = {}
    for i in participants:
        # Each participant generates a secret random nonce k_i
        nonces_k[i] = random.randint(1, q - 1)
        # And computes their public commitment R_i = g^k_i
        commitments_R[i] = pow(g, nonces_k[i], p)
        print(f"Party {i}: Generates secret nonce k_{i} and broadcasts commitment R_{i} = {commitments_R[i]}")
    
    # -- Round 2: Signature Share Generation --
    print("\n-- Round 2: Signature Share Generation --")
    # Each participant computes the group commitment R
    group_commitment_R = 1
    for i in participants:
        group_commitment_R = (group_commitment_R * commitments_R[i]) % p
    print(f"All participants compute the group commitment R = {group_commitment_R}")
    
    # Each participant computes the challenge e = H(R, Y, m)
    # In a real system, H is a cryptographic hash function like SHA-256.
    h = hashlib.sha256()
    h.update(str(group_commitment_R).encode())
    h.update(str(master_public_key_Y).encode())
    h.update(message.encode())
    challenge_e = int(h.hexdigest(), 16) % q
    print(f"All participants compute the challenge e = H(R, Y, m) = {challenge_e}")
    
    # Each participant computes their signature share s_i
    signature_shares_s = {}
    for i in participants:
        # Get the private share for this participant
        x_i = private_shares[i]
        # Get the nonce for this participant
        k_i = nonces_k[i]
        # Calculate the Lagrange coefficient for this participant
        lambda_i = get_lagrange_coeff(i, participants, q)
        # Compute the share
        s_i = (k_i + challenge_e * lambda_i * x_i) % q
        signature_shares_s[i] = s_i
        print(f"Party {i}: Computes signature share s_{i} = {s_i}")
    print("-" * 35 + "\n")

    # --- Phase 3: Aggregation and Verification ---
    print("--- Phase 3: Aggregation and Verification ---")
    
    # Aggregation: Sum the shares to get the final signature scalar s
    final_signature_s = sum(signature_shares_s.values()) % q
    
    print(f"Aggregated Signature: (R, s)")
    print(f"  R = {group_commitment_R}")
    print(f"  s = {final_signature_s}\n")
    
    # Verification
    print("-- Verification --")
    print("To verify, we check if g^s == R * Y^e (mod p)")
    
    # Left side of the equation: g^s
    lhs = pow(g, final_signature_s, p)
    
    # Right side of the equation: R * Y^e
    rhs = (group_commitment_R * pow(master_public_key_Y, challenge_e, p)) % p
    
    print("\nFinal Equation Numbers:")
    print(f"  g = {g}")
    print(f"  s = {final_signature_s}")
    print(f"  R = {group_commitment_R}")
    print(f"  Y = {master_public_key_Y}")
    print(f"  e = {challenge_e}\n")

    print(f"Left side (g^s mod p):  {lhs}")
    print(f"Right side (R*Y^e mod p): {rhs}")
    
    if lhs == rhs:
        print("\nResult: SUCCESS! The signature is valid.")
    else:
        print("\nResult: FAILURE! The signature is invalid.")

if __name__ == '__main__':
    main()
<<<SUCCESS! The signature is valid.>>>