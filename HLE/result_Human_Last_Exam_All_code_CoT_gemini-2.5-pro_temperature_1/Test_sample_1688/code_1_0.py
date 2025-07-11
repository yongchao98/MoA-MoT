import os

def xor_cipher(data, key):
    """A simple XOR cipher demonstration."""
    key_bytes = key.encode('utf-8')
    data_bytes = data.encode('utf-8')
    
    # Repeat the key to match the length of the data
    long_key = (key_bytes * (len(data_bytes) // len(key_bytes) + 1))[:len(data_bytes)]
    
    # Perform XOR operation byte by byte
    ciphered_bytes = bytes([b ^ k for b, k in zip(data_bytes, long_key)])
    
    # Try to decode back to a string
    try:
        return ciphered_bytes.decode('utf-8')
    except UnicodeDecodeError:
        return ciphered_bytes

# --- Scenario ---

# 1. A secret message and a secret key are established.
plaintext = "The launch code is 1234."
encryption_key = "SECRET_KEY_ALPHA"

print(f"Original Message: '{plaintext}'")
print(f"Encryption Key: '{encryption_key}'")
print("-" * 30)

# 2. The message is encrypted using the key.
# In a real system, only the sender and receiver have the key.
ciphertext_bytes = xor_cipher(plaintext, encryption_key).encode('utf-8')
print(f"Encrypted Ciphertext (as bytes): {ciphertext_bytes}")
print("-" * 30)


# 3. An adversary gains knowledge of BOTH the ciphertext AND the encryption key.
# This is the scenario described in the user's question.
print("Adversary's Turn:")
print("The adversary has intercepted the ciphertext and knows the key.")
adversary_known_ciphertext = ciphertext_bytes
adversary_known_key = encryption_key

print(f"Adversary knows Ciphertext: {adversary_known_ciphertext}")
print(f"Adversary knows Key: '{adversary_known_key}'")
print("\nAdversary's Decryption Process:")

# 4. The adversary uses the known key to decrypt the message.
# The decryption equation is: Plaintext = Ciphertext XOR Key
# The code will now perform this operation.

# To show the equation, we will print the values being used.
# Note: We convert back to string for the demo, the real operation is on bytes.
decrypted_text_by_adversary = xor_cipher(adversary_known_ciphertext.decode('utf-8', 'ignore'), adversary_known_key)

print(f"Equation: Plaintext = Ciphertext XOR Key")
print(f"Result: '{decrypted_text_by_adversary}' = {adversary_known_ciphertext} XOR '{adversary_known_key}'")
print("\nAdversary successfully recovered the message:")
print(f"Decrypted Message: '{decrypted_text_by_adversary}'")
print("\nConclusion: The system is compromised because the key was known.")
