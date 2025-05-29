def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_base = ord('A') if char.isupper() else ord('a')
            decrypted_char = chr((ord(char) - shift_base - shift) % 26 + shift_base)
            decrypted_text.append(decrypted_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "SOHDVH FKHFN WKH SURMHFW JXWHQEHUJ ZHE SDJHV IRU FXUUHQW GRQDWLRQ PHWKRGV DQG DGGUHVVHV."

for shift in range(1, 26):
    decrypted_text = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted_text}")