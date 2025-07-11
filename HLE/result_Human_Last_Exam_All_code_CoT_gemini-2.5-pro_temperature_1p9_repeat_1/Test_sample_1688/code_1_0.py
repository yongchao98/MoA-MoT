import os

def xor_cipher(text, key):
    """A simple XOR cipher function."""
    key_stream = (key[i % len(key)] for i in range(len(text)))
    return bytes([a ^ b for a, b in zip(text, key_stream)])

def main():
    """
    Demonstrates that knowledge of the key breaks the encryption.
    """
    # 1. A secret message is defined.
    plaintext = b"The meeting is at dawn."

    # 2. A secret key is generated. In a real system, this is kept secret.
    key = os.urandom(16)

    # 3. The message is encrypted using the key.
    #    Only the ciphertext is sent over an insecure channel.
    ciphertext = xor_cipher(plaintext, key)

    print(f"Original Message: {plaintext.decode()}")
    print(f"Encryption Key  : {key.hex()}")
    print(f"Transmitted Data: {ciphertext.hex()}")
    print("-" * 20)
    print("An adversary intercepts the data and knows the key.")
    print("-" * 20)

    # 4. The adversary's action:
    #    The adversary has the ciphertext and has also obtained the key.
    #    They use the same (publicly known) algorithm to decrypt.
    adversary_decrypted_text = xor_cipher(ciphertext, key)

    print("Adversary's Result:")
    print(f"Known Key       : {key.hex()}")
    print(f"Intercepted Data: {ciphertext.hex()}")
    print(f"Decrypted Message : {adversary_decrypted_text.decode()}")

    if adversary_decrypted_text == plaintext:
        print("\nConclusion: Security is compromised. The adversary successfully decrypted the message.")
    else:
        print("\nConclusion: Something went wrong, this shouldn't happen.")


if __name__ == "__main__":
    main()
