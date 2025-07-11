import ecdsa
import hashlib

def solve_crypto_puzzle():
    """
    Analyzes properties of digital signatures and demonstrates ECDSA malleability.
    """
    # 1. Generate an ECDSA key pair (using the SECP256k1 curve)
    # sk is the secret key, vk is the public key
    sk = ecdsa.SigningKey.generate(curve=ecdsa.SECP256k1)
    vk = sk.get_verifying_key()
    
    # 2. Define a message and sign it
    message = b"An example message"
    message_hash = hashlib.sha256(message).digest()
    
    # Generate the original signature. The result is a byte string.
    original_sig_bytes = sk.sign_digest(message_hash, sigencode=ecdsa.util.sigencode_der)
    
    # Decode the signature into its (r, s) integer components for inspection
    r, s = ecdsa.util.sigdecode_der(original_sig_bytes, vk.pubkey.order)
    
    print("--- Original Signature Details ---")
    print(f"Message: '{message.decode()}'")
    print(f"Public Key (uncompressed): {vk.to_string('uncompressed').hex()}")
    print(f"Original Signature (r): {hex(r)}")
    print(f"Original Signature (s): {hex(s)}")

    # 3. Verify the original signature to confirm it's valid
    try:
        is_valid = vk.verify_digest(original_sig_bytes, message_hash)
        assert is_valid
        print("\nVerification of original signature (r, s): SUCCESS")
    except ecdsa.BadSignatureError:
        print("\nVerification of original signature (r, s): FAILED")

    # 4. Create a new, different signature by manipulating 's'
    # The new s' will be equal to -s mod n, where n is the order of the curve's base point
    curve_order = vk.pubkey.order
    s_prime = curve_order - s
    
    # Encode the new signature (r, s') into the DER format
    malleated_sig_bytes = ecdsa.util.sigencode_der(r, s_prime, curve_order)
    
    print("\n--- Malleated Signature Details ---")
    print(f"Created a new signature by computing s' = n - s")
    print(f"Malleated Signature (r): {hex(r)}")
    print(f"Malleated Signature (s'): {hex(s_prime)}")

    # 5. Verify the new, malleated signature
    try:
        is_valid_malleated = vk.verify_digest(malleated_sig_bytes, message_hash)
        assert is_valid_malleated
        print("\nVerification of malleated signature (r, s'): SUCCESS")
    except ecdsa.BadSignatureError:
        print("\nVerification of malleated signature (r, s'): FAILED")

    print("\n\n--- Conclusion of Analysis ---")
    print("A. Create a new signature sig' for the same message with negligible probability? -> FALSE. We just created a valid new signature with probability 1.")
    print("B. Recover the secret key sk with negligible probability? -> TRUE. This is the fundamental security assumption (ECDLP hardness).")
    print("C. Forge a signature for a new message with non-negligible probability? -> FALSE. This is what a secure scheme like ECDSA prevents.")
    print("D. Recover message m from signature sig for all schemes with negligible probability? -> FALSE. Not true for schemes with message recovery.")

# Execute the demonstration and analysis
solve_crypto_puzzle()