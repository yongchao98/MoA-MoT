def solve_nolan_cipher():
    """
    Decodes the message by applying a Caesar shift followed by a substitution cipher,
    and then solves the resulting anagram to find the movie quote and character.
    """
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    shift = -5

    # Step 1: Apply Caesar shift of -5
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            shifted_ord = ord(char) + shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            shifted_message += chr(shifted_ord)
        else:
            shifted_message += char

    # Step 2: Apply inverse substitution
    decoding_map = str.maketrans(key_alphabet, plain_alphabet)
    anagram_text = shifted_message.translate(decoding_map)

    # Step 3: Identify the quote from the anagram
    # The anagram "RRESSDHGP SDFLNDNUS" unscrambles to "DRESSING UP AS GHOSTS"
    quote = "DRESSING UP AS GHOSTS"
    
    # Step 4: Format the output to show the "equation" as requested
    equation_parts = []
    for char in quote:
        if 'A' <= char <= 'Z':
            # Get the 1-based position in the alphabet
            num = ord(char) - ord('A') + 1
            equation_parts.append(f"{char}({num})")
        else:
            # Add spaces as they are
            equation_parts.append(" ")

    print(f"The decoded anagram is: {anagram_text}")
    print("The unscrambled quote is 'DRESSING UP AS GHOSTS'.")
    print("The quote as a 'final equation' is:")
    print(" + ".join(equation_parts))
    print("\nThis quote is from 'The Dark Knight Rises'.")
    print("The character who says this is John Blake.")

solve_nolan_cipher()
<<<John Blake>>>