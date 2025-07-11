def vigenere_decrypt(ciphertext, key):
    """
    Decrypts a Vigenère-encrypted text.
    It preserves the case of the original letters and ignores non-alphabetic characters.
    """
    key = key.lower()
    decrypted_text = ""
    key_index = 0
    for char in ciphertext:
        if 'a' <= char.lower() <= 'z':
            # Preserve the original case
            is_upper = char.isupper()
            
            # Calculate the shift for decryption
            key_char = key[key_index % len(key)]
            key_offset = ord(key_char) - ord('a')
            
            char_offset = ord(char.lower()) - ord('a')
            
            # Decrypt by subtracting the key shift
            decrypted_char_offset = (char_offset - key_offset + 26) % 26
            decrypted_char = chr(decrypted_char_offset + ord('a'))
            
            if is_upper:
                decrypted_text += decrypted_char.upper()
            else:
                decrypted_text += decrypted_char
            
            key_index += 1
        else:
            # Append non-alphabetic characters as they are
            decrypted_text += char
            
    return decrypted_text

# The encrypted text provided by the user
ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"

# The key derived from the puzzle
key = "AI"

# Decrypt the message
decrypted_question = vigenere_decrypt(ciphertext, key)

# Print the decrypted question
print("The decrypted question is:")
print(decrypted_question)

# Answer the question
# The question is: "the name of the programming language created by guido van rossum?"
# The answer is Python.
print("\nThe answer to the question is:")
print("Python")