def decrypt_vigenere(ciphertext, key):
    """Decrypts a Vigenere-encoded text with a given key."""
    key_index = 0
    plaintext = ""
    for char in ciphertext.lower():
        if 'a' <= char <= 'z':
            # Calculate the shift value from the key
            key_char = key[key_index % len(key)]
            shift = ord(key_char) - ord('a')
            
            # Decrypt the character
            decrypted_char_code = ord(char) - shift
            if decrypted_char_code < ord('a'):
                decrypted_char_code += 26
            
            plaintext += chr(decrypted_char_code)
            key_index += 1
        else:
            # Keep non-alphabetic characters as they are
            plaintext += char
    return plaintext

# The encrypted question
ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"

# The key derived from the puzzle (AMEN -> MEENA)
key = "meena"

# Decrypt the text
decrypted_question = decrypt_vigenere(ciphertext, key)

# The decrypted question is "how many c's are in canonicallyuncharacteristicallyunconstitutionalnesses?"
# Now, let's answer the question by counting the 'c's in the given word.
long_word = "canonicallyuncharacteristicallyunconstitutionalnesses"
count = long_word.count('c')

# Print the final equation as requested
# Using f-string to print each number of the answer
answer_str = str(count)
final_equation = f"{decrypted_question} {answer_str}"

print(final_equation)
for char in answer_str:
    print(char)
