import os

def generate_otp_key(length):
    """Generates a secure, random key of a specific length."""
    return os.urandom(length)

def otp_encrypt(plaintext_bytes, key):
    """Encrypts plaintext using OTP by XORing with the key."""
    return bytes([p ^ k for p, k in zip(plaintext_bytes, key)])

def otp_decrypt(ciphertext_bytes, key):
    """Decrypts ciphertext using OTP by XORing with the same key."""
    # Encryption and decryption are the same operation in OTP
    return bytes([c ^ k for c, k in zip(ciphertext_bytes, key)])

def main():
    """Demonstrates the OTP process for a DNS query."""
    # A sample DNS query we want to protect
    dns_query = "www.example.com"
    plaintext_bytes = dns_query.encode('utf-8')

    # 1. Generate a one-time pad (key) of the same length as the message.
    # This key must be truly random, used only once, and shared securely.
    key = generate_otp_key(len(plaintext_bytes))

    # 2. Encrypt the DNS query.
    ciphertext_bytes = otp_encrypt(plaintext_bytes, key)

    # 3. Decrypt the message using the same key.
    decrypted_bytes = otp_decrypt(ciphertext_bytes, key)
    decrypted_text = decrypted_bytes.decode('utf-8')

    print("--- One-Time Pad (OTP) Demonstration ---")
    print(f"Original DNS Query: {dns_query}\n")

    # The "equation" is: Plaintext XOR Key = Ciphertext
    # We will print the numerical byte values for each part.
    print("Plaintext bytes (P):")
    print(' '.join(map(str, plaintext_bytes)))
    print("\nOne-Time Pad key (K):")
    print(' '.join(map(str, key)))
    print("\nCiphertext bytes (C = P ^ K):")
    print(' '.join(map(str, ciphertext_bytes)))

    print("\n--- Decryption ---")
    print(f"Decrypted Text: {decrypted_text}")

    # Verification
    if decrypted_bytes == plaintext_bytes:
        print("\nVerification successful: Decrypted text matches original plaintext.")
    else:
        print("\nVerification failed.")

if __name__ == "__main__":
    main()