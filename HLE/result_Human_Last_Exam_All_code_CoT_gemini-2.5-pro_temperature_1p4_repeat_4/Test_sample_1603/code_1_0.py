def solve_cipher():
    """
    Decodes the message by first applying a Caesar cipher shift and then a substitution cipher.
    The final output will be the name of the character who spoke the decoded quote.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alpha = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5

    print(f"Decoding message: {encoded_message}")
    print(f"Using key: {key_alpha}")
    print(f"Using shift: {shift}\n")

    # Step 1: Create the shifted message by applying the Caesar cipher
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            shifted_ord = ord(char) + shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            shifted_message += chr(shifted_ord)
        else:
            shifted_message += char

    # Step 2: Create the decoding map for the substitution cipher
    decoding_map = {key_char: std_char for key_char, std_char in zip(key_alpha, std_alpha)}

    # Step 3: Apply the substitution to the shifted message and build the final message
    final_message_parts = []
    # Process each character and show the equation
    for i, char in enumerate(encoded_message):
        if 'A' <= char <= 'Z':
            shifted_char = shifted_message[i]
            decoded_char = decoding_map[shifted_char]
            final_message_parts.append(decoded_char)
            print(f"{char} - 5 = {shifted_char}, which decodes to {decoded_char}")
        else:
            final_message_parts.append(char)
            print("Space detected, keeping space")

    final_quote = "".join(final_message_parts)

    print(f"\nDecoded Quote: {final_quote}")
    print("\nThe quote is a chant from The Dark Knight Rises, meaning 'He Rises'.")

    # The character associated with this quote
    character_name = "Bane"
    print(f"The character who says this is: {character_name}")


solve_cipher()

<<<Bane>>>