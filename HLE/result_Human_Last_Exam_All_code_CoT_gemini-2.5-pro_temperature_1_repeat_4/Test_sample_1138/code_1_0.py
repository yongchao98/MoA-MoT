import hashlib

# This is a TOY digital signature scheme, created for demonstration purposes only.
# It is completely insecure and should not be used in practice.

class InsecureSignatureScheme:
    """
    This class simulates a digital signature scheme that is universally forgeable.
    The vulnerability is that the secret key can be trivially derived from the public key.
    """

    def keygen(self):
        """Generates a public/secret key pair."""
        # In a real scheme, sk would be a large random number.
        sk = 12345
        # The public key is insecurely derived, making the secret key discoverable.
        pk = sk * 2
        return pk, sk

    def sign(self, sk: int, message: str) -> str:
        """Signs a message with the secret key."""
        message_with_key = f"{message}:{sk}"
        # Use SHA256 as the hash function.
        signature = hashlib.sha256(message_with_key.encode()).hexdigest()
        return signature

    def verify(self, pk: int, message: str, signature: str) -> bool:
        """Verifies a signature against the public key and message."""
        # This verification process reveals the vulnerability.
        # Anyone can derive the 'secret' key from the public key.
        derived_sk = pk // 2
        
        # Re-create the signature to verify it.
        expected_signature = self.sign(derived_sk, message)
        
        return signature == expected_signature

def adversary_forge(pk: int, message: str) -> str:
    """
    An adversary who, given only the public key and a target message,
    can forge a valid signature.
    """
    print(f"Adversary attempts to forge a signature for message: '{message}'")
    
    # The adversary exploits the vulnerability to derive the secret key.
    derived_sk = pk // 2
    print(f"Adversary derived secret key '{derived_sk}' from public key '{pk}'.")
    
    # The adversary can now sign any message.
    forged_signature = InsecureSignatureScheme().sign(derived_sk, message)
    print(f"Adversary created forged signature: {forged_signature}")
    
    return forged_signature

# --- Demonstration ---

# 1. Setup the scheme and generate keys.
scheme = InsecureSignatureScheme()
public_key, secret_key = scheme.keygen()
print(f"System generated Public Key: {public_key}, Secret Key: {secret_key}\n")

# 2. Define a target message for the adversary.
target_message = "Transfer 1 million dollars to Eve"

# 3. The adversary performs the universal forgery.
#    They only need the public_key and the message they want to sign.
forged_signature = adversary_forge(public_key, target_message)
print("\n--- Verification ---")

# 4. Anyone can now verify the forged signature, and it will be accepted as valid.
is_valid = scheme.verify(public_key, target_message, forged_signature)

print(f"Verification of forged signature for message '{target_message}': {is_valid}")

if is_valid:
    print("\nConclusion: The forgery was successful. This demonstrates a scheme that is 'universally forgeable',")
    print("which fits the description in option C, as such a scheme is a type of 'existentially forgeable' scheme.")
