def vigenere_decrypt(ciphertext, key):
    """Decrypts ciphertext using a Vigenere cipher."""
    key_len = len(key)
    plaintext = ""
    key_index = 0
    
    for char in ciphertext:
        if 'a' <= char <= 'z':
            key_char = key[key_index % key_len].lower()
            key_shift = ord(key_char) - ord('a')
            
            cipher_val = ord(char) - ord('a')
            plain_val = (cipher_val - key_shift + 26) % 26
            
            plaintext += chr(plain_val + ord('a'))
            key_index += 1
        else:
            # Keep non-alphabetic characters as they are
            plaintext += char
            
    return plaintext

# The ciphertext provided by the user.
ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"

# The key found through cryptanalysis, as the puzzle was a likely red herring.
key = "MEEAMEEMENAMEMEENAGLEEYES"

# Decrypt the message to reveal the question.
decrypted_question = vigenere_decrypt(ciphertext, key)
print(f"Decrypted Question: {decrypted_question}\n")

# Now, we answer the decrypted question.
# The question is: "how many i's are in the international phonetic alphabet representation of yes?"

# 1. Find the IPA for "yes".
ipa_of_yes = "/jɛs/"

# 2. Count the number of 'i's in the IPA string.
# The characters are 'j', 'ɛ' (epsilon), and 's'. None of them is 'i'.
count_of_i = ipa_of_yes.count('i')

# 3. Print the final answer and the reasoning.
print("To answer the question:")
print(f"1. The International Phonetic Alphabet (IPA) representation of 'yes' is '{ipa_of_yes}'.")
print(f"2. We need to count the number of times the letter 'i' appears in '{ipa_of_yes}'.")
print(f"3. The count is {count_of_i}.")
print(f"\nFinal Answer: {count_of_i}")

# The format for the final answer.
# The final equation is simply counting the occurrences of 'i' in the string "/jɛs/".
print(f"The number of 'i's in '{ipa_of_yes}' is {count_of_i}.")
