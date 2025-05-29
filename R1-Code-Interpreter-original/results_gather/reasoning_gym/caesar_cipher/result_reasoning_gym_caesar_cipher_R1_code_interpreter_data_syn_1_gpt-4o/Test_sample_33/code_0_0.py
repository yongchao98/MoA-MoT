def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
            elif char.islower():
                if shifted < ord('a'):
                    shifted += 26
            decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

encrypted_text = "DRO CKSN WB"
for shift in range(1, 26):
    decrypted_text = decrypt_caesar_cipher(encrypted_text, shift)
    print(f"Shift {shift}: {decrypted_text}")