import ecdsa
import hashlib

# This script demonstrates the malleability of ECDSA, which proves that
# statement A in the question is false.

# 1. Generate an ECDSA key pair (using the SECP256k1 curve)
sk = ecdsa.SigningKey.generate(curve=ecdsa.SECP256k1)
pk = sk.get_verifying_key()
curve_order = sk.curve.order

# 2. Define a message to be signed
message = b"An existentially forgeable signature scheme."
message_hash = hashlib.sha256(message).digest()

# 3. Sign the message to get an original signature (r, s)
# The library's sign() method returns the signature in DER format. We decode it.
signature_der = sk.sign(message_hash, hashfunc=hashlib.sha256)
r, s = ecdsa.util.sigdecode_string(signature_der, curve_order)

print("--- Original Signature (sig) ---")
print(f"Message: {message.decode()}")
print(f"r = {r}")
print(f"s = {s}")

# 4. Verify the original signature is valid
try:
    pk.verify(signature_der, message_hash, hashfunc=hashlib.sha256)
    print("\nVerification of original signature (sig): SUCCESS")
except ecdsa.BadSignatureError:
    print("\nVerification of original signature (sig): FAILED")

# 5. Create a new, malleated signature (sig')
# The new signature has the same 'r' but a new s' = n - s, where n is the curve order.
s_malleated = curve_order - s

# Re-encode the malleated signature to DER format for verification
malleated_sig_der = ecdsa.util.sigencode_string(r, s_malleated, curve_order)


print("\n--- Malleated Signature (sig') ---")
print("A new signature sig' = (r, s') is created from sig = (r, s).")
print(f"The equation is: s' = n - s")
print(f"n (curve order) = {curve_order}")
print(f"s (original)    = {s}")
print(f"s' (malleated)  = {curve_order} - {s} = {s_malleated}")

# 6. Verify the new signature is also valid for the same message
try:
    pk.verify(malleated_sig_der, message_hash, hashfunc=hashlib.sha256)
    print("\nVerification of malleated signature (sig'): SUCCESS")
except ecdsa.BadSignatureError:
    print("\nVerification of malleated signature (sig'): FAILED")

print("\n--- Conclusion from Code ---")
print("Both the original signature and the new, malleated signature are valid for the same message.")
print("This shows an adversary can create a new signature with probability 1, disproving statement A.")