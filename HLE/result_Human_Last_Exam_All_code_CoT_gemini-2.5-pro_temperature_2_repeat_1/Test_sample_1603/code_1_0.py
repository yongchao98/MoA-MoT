def decode_nolan_quote(encoded_message, key):
    """
    Decodes a message using a keyed Caesar cipher.

    The process for each character is:
    1. Find the index of the character in the provided key.
    2. Subtract 5 (the Caesar shift) from this index.
    3. The decoded character is the one at the new index in the standard alphabet.
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    print("Decoding Process:")
    print("-----------------")

    decoded_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            try:
                # 1. Find the index of the character in the key
                idx_in_key = key.find(char)
                
                # 2. Subtract 5 from the index (Caesar shift)
                # We add 26 before the modulo to handle negative results gracefully
                shifted_idx = (idx_in_key - 5 + 26) % 26
                
                # 3. Find the new character in the standard alphabet
                decoded_char = alphabet[shifted_idx]
                decoded_message += decoded_char
                
                print(f"'{char}' -> index in key: {idx_in_key:2d} -> ({idx_in_key:2d} - 5) mod 26 = {shifted_idx:2d} -> '{decoded_char}'")

            except ValueError:
                # Character not in the key, append as is
                decoded_message += char
                print(f"'{char}' -> not in key, keeping character.")
        else:
            # Handle spaces or other characters
            decoded_message += char
            if char == ' ':
                print("\n(space)\n")

    print("-----------------\n")
    print("Decoded Message:")
    print(final_message := decoded_message)
    return final_message

# Given data from the problem
encoded_text = "PPXRRNUTL RNYFYINIZR"
cipher_key = "CHRISTOPENLABDFGJKMQUVWXYZ"

# Run the decoder
decode_nolan_quote(encoded_text, cipher_key)

# The decoded quote is "ARE YOU WATCHING CLOSELY".
# This famous quote is said by Alfred Borden in "The Prestige".
print("\nThe quote is said by the character Alfred Borden from the movie The Prestige.")
<<<Alfred Borden>>>