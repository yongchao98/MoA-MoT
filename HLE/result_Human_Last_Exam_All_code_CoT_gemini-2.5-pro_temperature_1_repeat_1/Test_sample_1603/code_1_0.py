import collections

def solve_cipher():
    """
    Solves a multi-step cipher to reveal a quote from a Christopher Nolan movie.
    The cipher involves an inverse Caesar shift followed by an inverse substitution.
    """
    # 1. Define the cipher parameters
    encoded_message = "PPXRRNUTL RNYFYINIZR"
    plain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # Construct the substitution alphabet from the keyword
    keyword = "CHRISTOPHERNOLAN"
    key_alphabet_part = ""
    for char in keyword.upper():
        if char not in key_alphabet_part and 'A' <= char <= 'Z':
            key_alphabet_part += char
    
    key_alphabet = key_alphabet_part
    for char in plain_alphabet:
        if char not in key_alphabet:
            key_alphabet += char
            
    shift = 5

    # 2. Step 1 of Decryption: Apply inverse Caesar shift (-5)
    shifted_message = ""
    for char in encoded_message:
        if 'A' <= char <= 'Z':
            # Shift character code back by 5, wrapping around the alphabet
            shifted_code = ord(char) - shift
            if shifted_code < ord('A'):
                shifted_code += 26
            shifted_message += chr(shifted_code)
        else:
            # Keep spaces and other characters as they are
            shifted_message += char

    # 3. Step 2 of Decryption: Apply inverse substitution
    # Create a map from the key_alphabet back to the plain_alphabet
    decryption_map = {key_alphabet[i]: plain_alphabet[i] for i in range(len(plain_alphabet))}
    
    final_message = ""
    for char in shifted_message:
        if char in decryption_map:
            final_message += decryption_map[char]
        else:
            final_message += char
            
    # 4. Print the final decoded message
    print(f"The decoded quote is: {final_message}")
    print("This is a quote from 'Inception'.")

solve_cipher()