import hashlib
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import ec
from cryptography.hazmat.primitives.asymmetric.utils import (
    decode_dss_signature, encode_dss_signature
)
from cryptography.exceptions import InvalidSignature

def demonstrate_ecdsa_malleability():
    """
    Demonstrates the malleability of ECDSA to evaluate the given statements.
    """
    # 1. Generate an ECDSA private key and the corresponding public key.
    # We use the SECP256R1 curve, a common standard.
    private_key = ec.generate_private_key(ec.SECP256R1())
    public_key = private_key.public_key()
    
    # The order 'n' of the curve's subgroup is needed for the calculation.
    curve_order_n = private_key.curve.order

    # 2. Define a message and sign it.
    message = b"This is the message to be signed."
    
    # ECDSA signing is typically performed on a hash of the message.
    digest = hashlib.sha256(message).digest()

    original_signature_der = private_key.sign(
        message,
        ec.ECDSA(hashes.SHA256())
    )

    # 3. Decode the signature from DER format to get the (r, s) integers.
    r, s = decode_dss_signature(original_signature_der)

    print("--- Original Signature Verification ---")
    print(f"Message: {message.decode()}")
    print(f"Signature component r: {r}")
    print(f"Signature component s: {s}")
    
    try:
        public_key.verify(
            original_signature_der,
            message,
            ec.ECDSA(hashes.SHA256())
        )
        print("Verification Successful: The original signature is valid.\n")
    except InvalidSignature:
        print("Verification Failed: The original signature is invalid.\n")
        return

    # 4. Demonstrate malleability.
    # A new, valid signature (r, s') can be created where s' = -s mod n.
    # In modular arithmetic, -s mod n is equivalent to n - s.
    s_prime = curve_order_n - s
    
    # This is the "equation" part of creating the new signature
    print("--- Malleable Signature Creation ---")
    print("A known property of ECDSA is that if (r, s) is a valid signature,")
    print("then (r, -s mod n) is also a valid signature for the same message.")
    print("We calculate s' = n - s (which is congruent to -s mod n).")
    print(f"Curve order n = {curve_order_n}")
    print(f"Original s    = {s}")
    print(f"New s' (n-s)  = {s_prime}")

    # Encode the new (r, s') pair back into DER format.
    malleable_signature_der = encode_dss_signature(r, s_prime)

    print("\n--- Malleable Signature Verification ---")
    print(f"New signature component r: {r}")
    print(f"New signature component s': {s_prime}")
    
    # 5. Verify the new, malleable signature. It should also be valid.
    try:
        public_key.verify(
            malleable_signature_der,
            message,
            ec.ECDSA(hashes.SHA256())
        )
        print("Verification Successful: The new malleable signature is also valid.\n")
    except InvalidSignature:
        print("Verification Failed: The new malleable signature is invalid.\n")

    # 6. Conclusion based on the demonstration.
    print("--- Conclusion ---")
    print("Statement A says a new signature can be created with 'no more than negligible probability'.")
    print("Our demonstration shows a new signature can be created with probability 1. Therefore, Statement A is FALSE.")
    print("\nStatement B says recovering the secret key has 'no more than negligible probability'.")
    print("This is the core security assumption of ECDSA, based on the hardness of the Elliptic Curve Discrete Logarithm Problem. This statement is TRUE.")

if __name__ == '__main__':
    demonstrate_ecdsa_malleability()