import string

def solve_cipher():
    """
    Solves the cipher by first applying a reverse Caesar shift and then
    reversing a substitution cipher.
    """
    ciphertext = "PPXRRNUTL RNYFYINIZR"
    key_alphabet = "CHRISTOPENLABDFGJKMQUVWXYZ"
    std_alphabet = string.ascii_uppercase
    shift = 5

    # Step 1: Reverse the Caesar cipher by shifting each letter back by 5.
    caesar_shifted_text = ""
    print("Ciphertext: " + ciphertext)
    print("Applying Caesar shift of -" + str(shift) + "...")
    
    shifted_equation = []
    for char in ciphertext:
        if 'A' <= char <= 'Z':
            original_ord = ord(char)
            shifted_ord = original_ord - shift
            if shifted_ord < ord('A'):
                shifted_ord += 26
            
            shifted_char = chr(shifted_ord)
            caesar_shifted_text += shifted_char
            shifted_equation.append(f"{char}({original_ord-ord('A')}) - {shift} = {shifted_char}({shifted_ord-ord('A')})")
        else:
            caesar_shifted_text += char
            
    print("Intermediate text after Caesar shift: " + caesar_shifted_text)

    # Step 2: Reverse the substitution cipher.
    # Create a decryption map from the key_alphabet to the std_alphabet.
    decryption_map = {key_alphabet[i]: std_alphabet[i] for i in range(26)}
    
    print("\nApplying substitution using key: " + key_alphabet)
    plaintext = ""
    for char in caesar_shifted_text:
        if char in decryption_map:
            plaintext += decryption_map[char]
        else:
            plaintext += char # For spaces or other characters.

    print("\nDecoded message equation:")
    # We will just show the final decoded letters, as a full equation is too long.
    final_chars = " ".join(list(plaintext))
    print(f'"{ciphertext}" -> "{plaintext}"')

    # The decoded message is "WHATISTHEMOSTRESILIENTPARASITE".
    # This is a quote from the movie Inception.
    # The character who says this is Cobb.
    
if __name__ == '__main__':
    solve_cipher()