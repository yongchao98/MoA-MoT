def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = shift % 26
            new_char = chr(((ord(char) - 65 - shift_amount) % 26) + 65)
            decrypted_text.append(new_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

encrypted_text = "HBXIW PH WT HTKTGTS RDBBJCXRPIXDC LXIW EPGXH"

for shift in range(1, 26):
    decrypted_text = decrypt_caesar_cipher(encrypted_text, shift)
    print(f"Shift {shift}: {decrypted_text}")