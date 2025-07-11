import string

def demonstrate_compromise(plaintext, key):
    """
    Demonstrates that knowing the encryption key breaks the security of a system.
    This function simulates a simple Caesar cipher.
    """
    print(f"Original Plaintext: '{plaintext}'")
    print(f"The Secret Key (shift value): {key}")
    print("-" * 30)

    # 1. Encryption by the legitimate user
    alphabet = string.ascii_lowercase
    shifted_alphabet = alphabet[key:] + alphabet[:key]
    translation_table = str.maketrans(alphabet, shifted_alphabet)
    ciphertext = plaintext.lower().translate(translation_table)

    print(f"Step 1: Message is encrypted into ciphertext.")
    print(f"Ciphertext: '{ciphertext}'")
    print("-" * 30)

    # 2. The Adversary's Action
    print("Step 2: An adversary obtains both the ciphertext AND the key.")
    print(f"Adversary has knowledge of:")
    print(f"  - The system (Caesar cipher)")
    print(f"  - The ciphertext ('{ciphertext}')")
    print(f"  - The key ({key})")
    print("-" * 30)

    # 3. Decryption by the Adversary (System is Compromised)
    # The adversary simply reverses the process.
    decryption_table = str.maketrans(shifted_alphabet, alphabet)
    revealed_plaintext = ciphertext.translate(decryption_table)

    print("Step 3: Adversary uses the key to decrypt the message.")
    print(f"Decrypted message: '{revealed_plaintext}'")
    print("-" * 30)

    # Conclusion
    if revealed_plaintext == plaintext:
        print("Conclusion: The system is COMPROMISED.")
        print("Because the key was known, the adversary successfully recovered the original message.")
        print("\nThis demonstrates that no system can remain secure if the encryption key is known.")
        print("The correct answer to the conceptual question is therefore F.")
    else:
        print("Something went wrong with the demonstration.")


# --- Simulation Parameters ---
# The original secret message.
original_message = "security depends on key secrecy"
# The secret key. If an adversary knows this, the system is broken.
secret_key = 13

# Run the demonstration
demonstrate_compromise(original_message, secret_key)
<<<F>>>