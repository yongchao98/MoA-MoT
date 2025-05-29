from collections import Counter

# The encrypted message
ciphertext = "KYV ZEMRCZUZKP FI LEVEWFITVRSZCZKP FW REP GIFMZJZFE FW KYZJ RXIVVDVEK JYRCC EFK MFZU KYV IVDRZEZEX GIFMZJZFEJ"

# Function to decrypt the text with a given shift
def decrypt_caesar(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
            else:
                if shifted < ord('a'):
                    shifted += 26
            decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

# Frequency analysis to find the most common letter
def frequency_analysis(text):
    # Remove spaces and count frequency of each letter
    text = text.replace(" ", "")
    frequency = Counter(text)
    return frequency

# Perform frequency analysis
freq = frequency_analysis(ciphertext)
most_common_letter = freq.most_common(1)[0][0]

# Assuming the most common letter in the ciphertext is 'E' in plaintext
# Calculate the shift
shift = ord(most_common_letter) - ord('E')

# Decrypt the message with the calculated shift
decrypted_message = decrypt_caesar(ciphertext, shift)
print(decrypted_message)