import string

def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char in string.ascii_uppercase:
            index = (ord(char) - ord('A') - shift) % 26
            decrypted_text.append(chr(ord('A') + index))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "ZW PFLI JKFDRTY ZJ FLK FW ZK DLJK SV DVEUVU"

for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")