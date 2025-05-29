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

# Try all possible shifts and print the results
for shift in range(26):
    decrypted_message = decrypt_caesar(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_message}")