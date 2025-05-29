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

# The encrypted message with shift 4
shift_4_text = "GUR VAINYVQVGL BE HARASBEPRNOVYVGL BS NAL CEBIVFVBA BS GUVF NTERRZRAG FUNYY ABG IBVQ GUR ERZNVAVAT CEBIVFVBAF"

# Decrypt using ROT13 (shift of 13)
decrypted_message = decrypt_caesar(shift_4_text, 13)
print(decrypted_message)