def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_char = chr((ord(char) - start + shift_amount) % 26 + start)
            decrypted_text.append(decrypted_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "UGRF RFCGP FCYBQ ZCLR BMUL MTCP RFCGP CJCARPGA RFGPRW QAGCLRGDGA KCL UCPC YZQMPZCB GL RPYLQACLBCLRYJ AYJASJYRGMLQ."
decrypted_text = caesar_decrypt(ciphertext, 2)
print(decrypted_text)