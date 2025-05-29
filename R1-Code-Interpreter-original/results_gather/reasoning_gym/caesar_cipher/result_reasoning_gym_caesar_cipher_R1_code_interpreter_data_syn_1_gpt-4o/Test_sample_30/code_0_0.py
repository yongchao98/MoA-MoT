def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
            elif char.islower():
                if shifted < ord('a'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "IUVEXOMNZ RGC OT IXKGZOTM ZNK VXUPKIZ IURRKIZOUT"

for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")