def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
            decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "KDIRGOSY KBKXE UTK VXKYKTZ"
for shift in range(1, 26):
    print(f"Shift {shift}: {decrypt_caesar_cipher(cipher_text, shift)}")