def solve_nolan_cipher():
    """
    Decodes the message by applying a Vigenere cipher with a custom alphabet
    and a specific key derived from the hints.
    """
    # 1. Define the elements of the cipher from the problem description.
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    custom_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    # The number 5 hints at the 5-letter key "NOLAN".
    vigenere_key = "NOLAN"

    # Create a mapping from each character to its index for efficient lookup.
    alpha_to_idx = {char: i for i, char in enumerate(custom_alphabet)}

    decoded_message = ""
    key_index = 0

    print(f"Decoding message: '{encoded_message}'")
    print(f"Using Vigenere Key: '{vigenere_key}'")
    print(f"Using Custom Alphabet: '{custom_alphabet}'")
    print("-" * 20)
    print("Decoding steps (Equation: Plaintext_Char = (Cipher_Char - Key_Char)):\n")

    # 2. Iterate through each character of the encoded message.
    # The message has a space, which should be preserved. We'll handle the words separately.
    words = encoded_message.split(" ")
    full_decoded_message = []

    for word in words:
        decoded_word = ""
        for char in word:
            # Get the index of the ciphertext character.
            cipher_char_index = alpha_to_idx[char]

            # Get the current key character and its index.
            key_char = vigenere_key[key_index % len(vigenere_key)]
            key_char_index = alpha_to_idx[key_char]

            # 3. Apply the Vigenere decoding formula.
            # Plaintext_Index = (Ciphertext_Index - Key_Index + 26) % 26
            plain_char_index = (cipher_char_index - key_char_index + 26) % 26

            # Find the resulting plaintext character from the custom alphabet.
            plain_char = custom_alphabet[plain_char_index]
            decoded_word += plain_char

            # Print the equation for the current step.
            print(f"'{plain_char}' = ('{char}' - '{key_char}') -> {plain_char_index} = ({cipher_char_index} - {key_char_index}) mod 26")

            # Move to the next character in the Vigenere key.
            key_index += 1
        full_decoded_message.append(decoded_word)

    final_quote = " ".join(full_decoded_message)

    print("-" * 20)
    print(f"Final Decoded Quote: {final_quote}")
    print("\nThe quote is from the movie 'Inception'.")
    print("The character who says this quote is Mal Cobb.")

solve_nolan_cipher()
<<<Mal>>>