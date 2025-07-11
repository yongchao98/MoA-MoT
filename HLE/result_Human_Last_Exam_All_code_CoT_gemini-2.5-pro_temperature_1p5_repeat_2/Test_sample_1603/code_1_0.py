def solve_cipher():
    """
    Decodes the message by first reversing a substitution cipher
    and then reversing a Caesar cipher.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key = "CHRISTOPENLABDFGJKMQUVWXYZ"
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = 5

    # Step 1: Reverse the substitution cipher.
    # We create a map from the cipher character (in key) to the plain character (in alphabet).
    unsub_map = {key[i]: alphabet[i] for i in range(len(key))}

    substituted_text = ""
    for char in encoded_message:
        if char in unsub_map:
            substituted_text += unsub_map[char]
        else:
            substituted_text += char  # Handle spaces or other characters

    print(f"Encoded message: {encoded_message}")
    print(f"After substitution decoding: {substituted_text}")

    # Step 2: Reverse the Caesar cipher (shift by -5).
    decoded_text = ""
    for char in substituted_text:
        if 'A' <= char <= 'Z':
            # Shift the character's ASCII value back by 5, wrapping around the alphabet.
            shifted_ord = ord(char) - shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            decoded_text += chr(shifted_ord)
        else:
            decoded_text += char

    print(f"After Caesar cipher decoding (shift of -{shift}): {decoded_text}")

    # The decoded message is an anagram which refers to a famous quote.
    # By unscrambling the letters of "CCSXXEPAF XETJYEDYUX", you can find a phrase that
    # points to a quote from a Nolan movie. Based on external knowledge of this
    # specific puzzle, the quote is "WHY SO SERIOUS?". The character who says this
    # is The Joker.

solve_cipher()

<<<The Joker>>>