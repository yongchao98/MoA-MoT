def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = shift % 26
            shifted = ord(char) - shift_amount
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

cipher_text = "DOLYL OHCL FVB HZRLK TY."
for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")