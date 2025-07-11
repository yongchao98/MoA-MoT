import base64

def encrypt(plaintext, key):
    """
    A simple XOR cipher. The key is repeated to match the length of the plaintext.
    This demonstrates a basic symmetric encryption algorithm.
    """
    key_bytes = key.encode('utf-8')
    plaintext_bytes = plaintext.encode('utf-8')
    
    encrypted_bytes = bytearray()
    for i in range(len(plaintext_bytes)):
        encrypted_bytes.append(plaintext_bytes[i] ^ key_bytes[i % len(key_bytes)])
        
    # Using base64 to make the output printable text
    return base64.b64encode(encrypted_bytes).decode('utf-8')

def decrypt(ciphertext, key):
    """
    The corresponding decryption function for the simple XOR cipher.
    Note: It's the same operation as encryption in this case.
    """
    key_bytes = key.encode('utf-8')
    ciphertext_bytes = base64.b64decode(ciphertext.encode('utf-8'))
    
    decrypted_bytes = bytearray()
    for i in range(len(ciphertext_bytes)):
        decrypted_bytes.append(ciphertext_bytes[i] ^ key_bytes[i % len(key_bytes)])
        
    return decrypted_bytes.decode('utf-8')

# --- Scenario ---
# The cryptographic system (encrypt/decrypt functions) is public.
# The encryption key MUST remain secret for the system to be secure.

# Let's define the secret elements
secret_message = "True security is impossible if the key is public."
secret_key = "SuperSecretKey123"

print(f"The original secret message is: '{secret_message}'")
print(f"The secret key is: '{secret_key}'")

# Encrypting the message
encrypted_message = encrypt(secret_message, secret_key)
print(f"The encrypted message is: '{encrypted_message}'")
print("-" * 20)

# --- Adversary's Perspective ---
# The question states the adversary has complete knowledge, including the key.
adversary_knows_protocol = True
adversary_knows_key = True
adversary_intercepted_message = encrypted_message

print("An adversary operates under the following conditions:")
print(f"  - Knows the protocol (the 'decrypt' function): {adversary_knows_protocol}")
print(f"  - Knows the key ('{secret_key}'): {adversary_knows_key}")
print(f"  - Intercepted the ciphertext: '{adversary_intercepted_message}'")
print("-" * 20)

# The adversary now uses their knowledge to break the encryption
if adversary_knows_protocol and adversary_knows_key:
    print("Adversary is attempting to decrypt the message...")
    # The adversary simply runs the public decryption algorithm with the known key
    adversary_decrypted_message = decrypt(adversary_intercepted_message, secret_key)
    
    print(f"\nMessage decrypted by adversary: '{adversary_decrypted_message}'")
    
    if adversary_decrypted_message == secret_message:
        print("\nConclusion: Security is compromised. The system cannot remain secure under these conditions.")
