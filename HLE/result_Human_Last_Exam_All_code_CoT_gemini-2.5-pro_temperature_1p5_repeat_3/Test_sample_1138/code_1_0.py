import hashlib
from ecdsa import SigningKey, NIST256p
from ecdsa.util import sigencode_string, sigdecode_string

def demonstrate_ecdsa_properties():
    """
    This function demonstrates properties of ECDSA to analyze the provided options.
    It shows that creating a new signature for a pre-existing message is possible
    (making statement A false), while reinforcing that this does not imply the secret key
    is recoverable (making statement B true).
    """
    # 1. Generate a standard ECDSA key pair.
    # The secret key (sk) must be kept private.
    # The public key (pk) can be shared.
    sk = SigningKey.generate(curve=NIST256p)
    pk = sk.verifying_key
    curve_order = sk.curve.order

    print("Step 1: Generated an ECDSA key pair.")
    # print(f"Secret Key (d): {sk.privkey.secret_multiplier}") # Keep this secret in real life!
    print(f"Public Key (hex): {pk.to_string('hex')}")
    print("-" * 20)

    # 2. Sign a message with the secret key.
    message = b"An example message for signing"
    signature = sk.sign(message, hashfunc=hashlib.sha256)

    # The signature consists of two integers, r and s.
    r, s = sigdecode_string(signature, curve_order)
    print("Step 2: Signed a message with the secret key.")
    print(f"Message: '{message.decode()}'")
    print(f"Original Signature (r, s):")
    print(f"  r = {r}")
    print(f"  s = {s}")
    print("-" * 20)

    # 3. Verify the original signature with the public key.
    is_valid_original = pk.verify(signature, message, hashfunc=hashlib.sha256)
    print(f"Step 3: Verification of the original signature: {'SUCCESS' if is_valid_original else 'FAIL'}")
    print("-" * 20)
    
    # 4. Create a new, 'malleable' signature from the original.
    # The forgery is that for a signature (r, s), the signature (r, -s mod n) is also valid.
    # Here n is the order of the curve's base point.
    s_malleable = curve_order - s
    malleable_signature = sigencode_string(r, s_malleable, curve_order)

    print("Step 4: Created a new 'malleable' signature from the original.")
    print("This is done by calculating s' = n - s, where n is the curve order.")
    print(f"Malleable Signature (r, s'):")
    print(f"  r = {r}")
    print(f"  s' = {s_malleable}")
    print("-" * 20)

    # 5. Verify the malleable signature with the same public key.
    is_valid_malleable = pk.verify(malleable_signature, message, hashfunc=hashlib.sha256)
    print(f"Step 5: Verification of the malleable signature: {'SUCCESS' if is_valid_malleable else 'FAIL'}")
    print("-" * 20)

    # 6. Final Analysis
    print("Conclusion from the demonstration:")
    print("We created a new, different signature for the same message with 100% success.")
    print("This demonstrates that for ECDSA, creating a new signature 'sig\'' is possible with NON-NEGLIGIBLE probability.")
    print("\nTherefore, Statement A is FALSE.")
    print("\nThis malleability does NOT allow us to find the secret key 'sk'. Key recovery security is based on the hardness of the ECDLP, a different problem.")
    print("Therefore, Statement B, which claims an adversary can only recover the secret key with negligible probability, remains TRUE.")

if __name__ == '__main__':
    demonstrate_ecdsa_properties()