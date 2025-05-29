def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_char = chr((ord(char) - start - shift_amount) % 26 + start)
            decrypted_text.append(decrypted_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "UGRF RFCGP FCYBQ ZCLR BMUL MTCP RFCGP CJCARPGA RFGPRW QAGCLRGDGA KCL UCPC YZQMPZCB GL RPYLQACLBCLRYJ AYJASJYRGMLQ."

# Trying shift 2
decrypted = decrypt_caesar_cipher(cipher_text, 2)
print(decrypted)