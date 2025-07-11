import hashlib
import random

# A helper function for modular inverse, used for Lagrange coefficients
def mod_inverse(a, m):
    """Computes the modular inverse of a modulo m using the Extended Euclidean Algorithm."""
    a = a % m
    return pow(a, -2, m) # Using Fermat's Little Theorem, a^(m-2) = a^(-1) mod m for prime m

def print_equation(parts, result, modulus):
    """Helper function to print equations clearly."""
    equation = " + ".join(map(str, parts))
    print(f"    = {equation} mod {modulus}")
    print(f"    = {sum(parts)} mod {modulus}")
    print(f"    = {result}\n")

def main():
    """
    Simulates a 2-out-of-3 FROST-like Threshold Signature Scheme.
    NOTE: This is a simplified simulation for educational purposes. It uses modular
    arithmetic instead of elliptic curves and is NOT cryptographically secure.
    """
    print("--- 1. System Setup ---")
    # q: A large prime number for our finite field
    q = 104729
    # n: Total number of participants
    n = 3
    # t: Threshold of signers required
    t = 2
    # The message to be signed
    message = "Cryptography is fascinating"

    print(f"Prime Modulus (q): {q}")
    print(f"Total Participants (n): {n}")
    print(f"Threshold (t): {t}")
    print(f"Message (m): '{message}'\n")

    print("--- 2. Key Generation (Distributed Key Generation) ---")
    # In a real system, this is a complex interactive protocol (VSS).
    # Here, we simulate it with a trusted dealer.
    # We create a secret polynomial of degree t-1: f(z) = a_0 + a_1*z + ...
    # The main secret key is f(0) = a_0.
    poly_coeffs = [random.randint(1, q-1) for _ in range(t)]
    secret_key = poly_coeffs[0]

    # Each participant 'i' gets a secret share f(i)
    secret_shares = {}
    for i in range(1, n + 1):
        share = sum(poly_coeffs[j] * (i ** j) for j in range(t)) % q
        secret_shares[i] = share

    # The public key is derived from the main secret key (in real ECC: Y = secret_key * G)
    public_key = secret_key # Simplified for this simulation

    print(f"The secret polynomial is f(z) = {poly_coeffs[0]} + {poly_coeffs[1]}*z")
    print(f"Main Secret Key (f(0)): {secret_key}")
    print(f"Public Key: {public_key}")
    print("Secret Shares distributed:")
    for i, share in secret_shares.items():
        print(f"  Participant {i} has secret share: {share}")
    print("")

    print("--- 3. Signing Protocol ---")
    # Let's assume Participants 1 and 2 will sign the message.
    signing_participants = [1, 2]
    print(f"Signing participants: {signing_participants}\n")

    # Pre-computation: Each signer computes their Lagrange coefficient.
    # L_i = product of (j / (j - i)) for j in signing_participants where j != i
    print("Step A: Calculate Lagrange Coefficients for each signer")
    lagrange_coeffs = {}
    for i in signing_participants:
        num, den = 1, 1
        for j in signing_participants:
            if i != j:
                num = (num * j) % q
                den = (den * (j - i)) % q
        lagrange_coeffs[i] = (num * mod_inverse(den, q)) % q
        print(f"  Lagrange Coefficient for P{i} (L_{i}): {lagrange_coeffs[i]}")
    print("")

    # --- ROUND 1: Commitment ---
    print("--- Round 1: Commitment ---")
    nonces = {}
    commitments = {}
    for i in signing_participants:
        # Each participant generates two private nonces
        d_i, e_i = random.randint(1, q-1), random.randint(1, q-1)
        nonces[i] = (d_i, e_i)
        # In ECC, commitments are D_i = d_i*G, E_i = e_i*G. We use the numbers directly.
        commitments[i] = (d_i, e_i)
        print(f"P{i} generates nonces (d_{i}, e_{i}): ({d_i}, {e_i})")
        print(f"P{i} broadcasts commitments (D_{i}, E_{i}): ({commitments[i][0]}, {commitments[i][1]})")
    print("")

    # --- ROUND 2: Signature Share Generation ---
    print("--- Round 2: Signature Share Generation ---")
    # All participants have received the commitments from Round 1.

    # Step B: Compute Binding Factors. This prevents malleability attacks.
    # b_i = H(i, message, {D_j, E_j for all j})
    print("Step B: Compute Binding Factors (rho_i)")
    binding_factors = {}
    hasher = hashlib.sha256()
    hasher.update(message.encode())
    for i in signing_participants:
        hasher.update(str(commitments[i][0]).encode())
        hasher.update(str(commitments[i][1]).encode())

    base_hash = hasher.digest()
    for i in signing_participants:
        participant_hasher = hashlib.sha256()
        participant_hasher.update(base_hash)
        participant_hasher.update(str(i).encode())
        binding_factors[i] = int(participant_hasher.hexdigest(), 16) % q
        print(f"  Binding Factor for P{i} (b_{i}): {binding_factors[i]}")
    print("")

    # Step C: Compute Group Commitment (R).
    # In ECC: R = sum(D_i + b_i * E_i). We simulate this with modular arithmetic.
    print("Step C: Compute Group Commitment (R)")
    R_parts = []
    for i in signing_participants:
        d_i, e_i = commitments[i]
        b_i = binding_factors[i]
        R_parts.append((d_i + b_i * e_i) % q)
    R = sum(R_parts) % q
    print(f"  R = (D_1 + b_1*E_1) + (D_2 + b_2*E_2) mod q")
    print(f"  R = ({commitments[1][0]} + {binding_factors[1]}*{commitments[1][1]}) + ({commitments[2][0]} + {binding_factors[2]}*{commitments[2][1]}) mod {q}")
    print(f"  Group Commitment R = {R}\n")

    # Step D: Compute Challenge Hash (c).
    # c = H(Public Key, R, Message)
    print("Step D: Compute Challenge Hash (c)")
    challenge_hasher = hashlib.sha256()
    challenge_hasher.update(str(public_key).encode())
    challenge_hasher.update(str(R).encode())
    challenge_hasher.update(message.encode())
    c = int(challenge_hasher.hexdigest(), 16) % q
    print(f"  c = H({public_key}, {R}, '{message}')")
    print(f"  Challenge c = {c}\n")

    # Step E: Each participant computes their signature share.
    # z_i = d_i + e_i*b_i + c*L_i*s_i
    print("Step E: Compute individual Signature Shares (z_i)")
    sig_shares = {}
    for i in signing_participants:
        d_i, e_i = nonces[i]
        b_i = binding_factors[i]
        L_i = lagrange_coeffs[i]
        s_i = secret_shares[i]
        z_i = (d_i + (e_i * b_i) + (c * L_i * s_i)) % q
        sig_shares[i] = z_i
        print(f"  Share from P{i} (z_{i}) = d_{i} + e_{i}*b_{i} + c*L_{i}*s_{i} mod q")
        print(f"  z_{i} = {d_i} + {e_i}*{b_i} + {c}*{L_i}*{s_i} mod q")
        print(f"  z_{i} = {z_i}")
    print("")

    # --- 4. Signature Aggregation ---
    print("--- 4. Signature Aggregation ---")
    # The aggregator (can be one of the signers) collects t signature shares.
    z = sum(sig_shares.values()) % q
    print("The final signature value 'z' is the sum of the shares:")
    print(f"  z = z_1 + z_2 mod q")
    print_equation(list(sig_shares.values()), z, q)

    final_signature = (R, z)
    print(f"*** Final Assembled Signature (R, z) = {final_signature} ***\n")

    # --- 5. Signature Verification ---
    print("--- 5. Signature Verification ---")
    # To verify, a third party checks if z*G == R + c*Y
    # In our simulation, this is z == R + c*secret_key (mod q)
    lhs = z
    rhs = (R + c * secret_key) % q
    print(f"Checking if z == R + c * secret_key (mod q)")
    print(f"  Left side (z): {lhs}")
    print(f"  Right side ({R} + {c} * {secret_key}) mod {q}: {rhs}")

    if lhs == rhs:
        print("\n*** SUCCESS: The signature is valid! ***")
    else:
        print("\n*** FAILURE: The signature is invalid! ***")

if __name__ == '__main__':
    main()