import os

def demonstrate_compromise(plaintext_message):
    """
    Demonstrates that if an encryption key is known, the system is compromised.
    This function uses a simple XOR cipher, the principle behind the One-Time Pad.
    """
    print(f"Original Plaintext: '{plaintext_message}'\n")

    # 1. Convert the plaintext message to a list of byte values (numbers).
    plaintext_bytes = list(plaintext_message.encode('utf-8'))
    print(f"Plaintext as Bytes: {plaintext_bytes}\n")

    # 2. Generate a key. In a real system, this must be secret.
    # For this demonstration, we will assume the adversary obtains this key.
    key = list(os.urandom(len(plaintext_bytes)))
    print(f"Encryption Key as Bytes: {key}\n")

    # 3. Encrypt the plaintext by XORing it with the key.
    ciphertext_bytes = [p ^ k for p, k in zip(plaintext_bytes, key)]
    print(f"Resulting Ciphertext as Bytes: {ciphertext_bytes}\n")

    # 4. --- Adversary's Perspective ---
    # The adversary has captured the ciphertext AND knows the key.
    print("--------------------------------------------------")
    print("Adversary has captured ciphertext and knows the key.")
    print("--------------------------------------------------\n")
    print("Adversary begins decryption...\n")

    # 5. The adversary decrypts the message by XORing the ciphertext with the known key.
    decrypted_bytes = []
    print("Decryption Calculation (Ciphertext Byte ^ Key Byte = Plaintext Byte):")
    for i in range(len(ciphertext_bytes)):
        c_byte = ciphertext_bytes[i]
        k_byte = key[i]
        p_byte = c_byte ^ k_byte
        decrypted_bytes.append(p_byte)
        # Outputting each number in the final equation as requested.
        print(f"  {c_byte:3d} ^ {k_byte:3d} = {p_byte:3d}  (which is the character '{chr(p_byte)}')")

    # 6. Convert the decrypted bytes back into a string to reveal the message.
    recovered_plaintext = bytes(decrypted_bytes).decode('utf-8')
    print(f"\nRecovered Plaintext: '{recovered_plaintext}'")
    print("\nConclusion: Security is compromised. The system is broken.")


# Run the demonstration.
demonstrate_compromise("The key is the secret.")
