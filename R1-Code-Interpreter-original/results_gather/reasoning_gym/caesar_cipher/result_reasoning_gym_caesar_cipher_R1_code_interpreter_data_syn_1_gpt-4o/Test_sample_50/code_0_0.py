def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            decrypted_text.append(chr((ord(char) - shift_amount - shift) % 26 + shift_amount))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "VHSG SGDHQ GDZCR ADMS CNVM NUDQ SGDHQ DKDBSQHB SGHQSX RBHDMSHEHB LDM VDQD ZARNQBDC HM SQZMRBDMCDMSZK BZKBTKZSHNMR"

for shift in range(1, 26):
    decrypted = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted}")