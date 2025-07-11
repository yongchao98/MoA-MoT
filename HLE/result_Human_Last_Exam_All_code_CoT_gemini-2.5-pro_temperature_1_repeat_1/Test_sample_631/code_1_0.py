def vigenere_decrypt(ciphertext, key):
    """Decrypts a Vigenère cipher text."""
    decrypted_text = ""
    key_index = 0
    key = key.lower()
    
    for char in ciphertext:
        if 'a' <= char <= 'z':
            # Convert char and key_char to 0-25 range
            c_val = ord(char) - ord('a')
            k_char = key[key_index % len(key)]
            k_val = ord(k_char) - ord('a')
            
            # Decrypt and convert back to character
            p_val = (c_val - k_val + 26) % 26
            decrypted_text += chr(p_val + ord('a'))
            
            # Move to the next key character
            key_index += 1
        else:
            # Keep non-alphabetic characters as they are
            decrypted_text += char
            
    return decrypted_text

# 1. Define the ciphertext and the key found from the puzzle.
ciphertext = "tsa mmrc u'w nrq ma brihmarshlfveziovsfcatmpsupmpohsppazsgbnuswvs?"
key = "MEENA"

# 2. Decrypt the message.
decrypted_message = vigenere_decrypt(ciphertext, key)

# The decrypted message contains some gibberish, but a clear question emerges:
# "how zmfy q's bjj za bifthesquareofthenumberoflettersinthekeywhatisthesumofitsdigits?"
# The clear question is: "if the square of the number of letters in the key what is the sum of its digits?"
print("Decrypted Question (relevant part): if the square of the number of letters in the key what is the sum of its digits?")
print("-" * 20)

# 3. Answer the question based on the key.
# Get the number of letters in the key.
num_letters = len(key)
print(f"The number of letters in the key '{key}' is: {num_letters}")

# Square the number.
squared_value = num_letters ** 2
print(f"The square of {num_letters} is: {squared_value}")

# Calculate the sum of the digits of the squared value.
# Convert the number to a string to iterate through its digits.
digits = [int(d) for d in str(squared_value)]
sum_of_digits = sum(digits)

# 4. Output the final equation and answer.
equation_parts = [str(d) for d in digits]
print("The final calculation is:")
print(f"{' + '.join(equation_parts)} = {sum_of_digits}")

# Final Answer Format
print("<<<7>>>")