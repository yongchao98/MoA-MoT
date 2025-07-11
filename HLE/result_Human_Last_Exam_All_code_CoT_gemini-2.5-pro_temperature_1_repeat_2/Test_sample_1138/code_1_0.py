import ecdsa
import hashlib

def solve_signature_problem():
    """
    Analyzes the properties of ECDSA and other signature schemes
    to answer the multiple-choice question.
    """

    # --- Setup for Demonstration (Proving A is False) ---
    # We will demonstrate ECDSA's signature malleability.

    # 1. Generate a predictable key pair for demonstration purposes.
    # We use a common curve, SECP256k1, used by Bitcoin.
    sk = ecdsa.SigningKey.from_secret_exponent(123456789, curve=ecdsa.SECP256k1)
    pk = sk.get_verifying_key()
    
    # n is the order of the curve, a very important number in the equations.
    n = sk.curve.order

    # 2. Sign a message.
    message = b"An example message for signing."
    message_hash = hashlib.sha256(message).digest()
    
    # The standard ecdsa.py library produces a DER-encoded signature.
    original_signature_der = sk.sign_digest(message_hash, sigencode=ecdsa.util.sigencode_der)
    
    # To see the components (r, s), we decode the signature.
    r, s = ecdsa.util.sigdecode_der(original_signature_der, n)

    print("--- ECDSA Malleability Demonstration ---")
    print(f"Message: '{message.decode()}'")
    print(f"Curve Order (n): {n}")
    print(f"Original Signature (r, s): ({r}, {s})")

    # 3. Verify the original signature. It should be valid.
    try:
        is_valid_original = pk.verify_digest(original_signature_der, message_hash)
        print("Verification of original signature: VALID")
    except ecdsa.BadSignatureError:
        print("Verification of original signature: INVALID")

    # 4. Create the malleable signature, sig' = (r, n-s).
    s_malleable = n - s
    malleable_signature_der = ecdsa.util.sigencode_der(r, s_malleable, n)
    print("\nCreating a new signature using malleability...")
    print(f"Equation for new signature: sig' = (r, n-s)")
    print(f"Value of n-s: {s_malleable}")
    print(f"Malleable Signature (r, n-s): ({r}, {s_malleable})")
    
    # 5. Verify the malleable signature. It should also be valid.
    try:
        is_valid_malleable = pk.verify_digest(malleable_signature_der, message_hash)
        print("Verification of malleable signature: VALID")
    except ecdsa.BadSignatureError:
        print("Verification of malleable signature: INVALID")

    print(f"Signatures are different: {original_signature_der != malleable_signature_der}")
    print("-" * 20)

    # --- Analysis of Answer Choices ---
    print("\n--- Analysis of Answer Choices ---")
    print("Based on the demonstration and cryptographic principles:\n")

    print("A. For ECDSA: Given m, sig, pk, a computationally bounded adversary can create a new, different signature sig' ... with no more than negligible probability.")
    print("   - RESULT: FALSE. The code above demonstrates that creating a new, valid signature for the same message is trivial due to malleability. The probability is 1, not negligible.")
    print("   - Equation used: sig' = (r, n - s)\n")

    print("B. For ECDSA: Given m, sig, pk, a computationally bounded adversary can recover the secret key sk with no more than negligible probability.")
    print("   - RESULT: TRUE. This is the core security promise of ECDSA. The difficulty of recovering 'sk' is based on the Elliptic Curve Discrete Logarithm Problem (ECDLP), which is believed to be computationally hard. Signature malleability does not help in solving for 'sk'.\n")

    print("C. For some existentially forgeable digital signature schemes: Only given m, pk, a computationally bounded adversary can generate sig' that (sig', m) is verifiable against pk with non-negligible probability.")
    print("   - RESULT: TRUE. This describes universal forgery. An existentially forgeable scheme is not necessarily universally forgeable. However, a completely broken scheme (e.g., where the secret key is public) would be both existentially and universally forgeable. Since the statement only requires this to be true for 'some' schemes, the statement holds.\n")

    print("D. For all existentially forgeable digital signature schemes: Only given sig, pk, a computationally bounded adversary can figure out m with no more than negligible probability.")
    print("   - RESULT: FALSE. The 'for all' makes this statement incorrect. A counter-example is 'textbook RSA' without hashing. It is existentially forgeable, but given a signature 'sig', anyone with the public key (e, N) can recover the message 'm'.")
    print("   - Equation for counter-example: m = (sig ^ e) mod N\n")
    
    print("--- Conclusion ---")
    print("Both statements B and C are true. Therefore, the correct option is E.")

if __name__ == '__main__':
    solve_signature_problem()