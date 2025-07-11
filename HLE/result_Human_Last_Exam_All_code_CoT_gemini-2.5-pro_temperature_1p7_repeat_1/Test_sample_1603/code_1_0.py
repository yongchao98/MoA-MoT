def solve_cipher():
    """
    Solves a two-step cipher on an encoded message.
    1. Applies a Caesar cipher shift.
    2. Applies a substitution cipher using a custom key.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    standard_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5

    # Step 1: Apply the Caesar cipher shift (-5)
    shifted_message = ""
    for char in encoded_message:
        if char in standard_alphabet:
            original_index = standard_alphabet.index(char)
            new_index = (original_index + shift) % 26
            shifted_message += standard_alphabet[new_index]
        else:
            shifted_message += char

    # Step 2: Apply the substitution cipher
    # Create a map to decode from the key_alphabet to the standard_alphabet
    decoding_map = {key_char: std_char for std_char, key_char in zip(standard_alphabet, key_alphabet)}

    final_decoded_message = ""
    for char in shifted_message:
        if char in decoding_map:
            # Find the index of the character in the key alphabet
            key_index = key_alphabet.index(char)
            # The decoded character is the one at the same index in the standard alphabet
            decoded_char = standard_alphabet[key_index]
            final_decoded_message += decoded_char
        else:
            final_decoded_message += char
            
    # Print the final result including the "equation" for each letter
    print(f"Decoded quote: {final_decoded_message}")
    print("\n--- Calculation Steps ---")
    
    # Isolate words for clarity
    words = final_decoded_message.split()
    
    # Print the equation for each letter of the decoded message
    for i, char in enumerate(final_decoded_message):
        if char.isalpha():
            # Find the original letter by reversing the process
            original_shifted_char = key_alphabet[standard_alphabet.index(char)]
            original_encoded_char = standard_alphabet[(standard_alphabet.index(original_shifted_char) - shift) % 26]

            # Find the numerical indexes
            original_encoded_index = standard_alphabet.index(original_encoded_char)
            shifted_index = (original_encoded_index + shift) % 26
            final_decoded_index = standard_alphabet.index(char)
            
            print(f"'{original_encoded_char}' (index {original_encoded_index}) "
                  f"-> Caesar shift by {shift} -> '{original_shifted_char}' (index {shifted_index}) "
                  f"-> Substitute via key -> '{char}' (index {final_decoded_index})")
        # Handle the space
        elif char.isspace():
            print("\n")


solve_cipher()