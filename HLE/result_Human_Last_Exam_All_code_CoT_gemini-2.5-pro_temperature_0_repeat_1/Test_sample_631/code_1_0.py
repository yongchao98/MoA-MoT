def vigenere_decrypt(ciphertext, key):
    """Decrypts a Vigenère-encrypted text."""
    decrypted_text = ""
    key_index = 0
    key = key.lower()
    for char in ciphertext:
        if 'a' <= char <= 'z':
            # Calculate the shift for decryption
            key_char = key[key_index % len(key)]
            key_shift = ord(key_char) - ord('a')
            
            cipher_val = ord(char) - ord('a')
            
            # Decrypt and handle wrap-around
            plain_val = (cipher_val - key_shift + 26) % 26
            
            decrypted_text += chr(plain_val + ord('a'))
            key_index += 1
        else:
            # Keep non-alphabetic characters as they are
            decrypted_text += char
    return decrypted_text

# The encrypted text provided by the user
ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"

# The key derived from the puzzle ("AMEN" -> "MEENA")
key = "meena"

# Decrypt the message to find the question
decrypted_question = vigenere_decrypt(ciphertext, key)
print(f"Decrypted Question: {decrypted_question}")

# The word of affirmation used to solve the puzzle for the key
word_of_affirmation = "AMEN"
print(f"The word of affirmation used to find the key was: {word_of_affirmation}")

# Answer the decrypted question
# The question asks for the number of 'i's in the word of affirmation.
answer = word_of_affirmation.lower().count('i')

# The final answer is the count of 'i's in "AMEN"
print(f"The number of 'i's in \"{word_of_affirmation}\" is {answer}.")

print("<<<0>>>")