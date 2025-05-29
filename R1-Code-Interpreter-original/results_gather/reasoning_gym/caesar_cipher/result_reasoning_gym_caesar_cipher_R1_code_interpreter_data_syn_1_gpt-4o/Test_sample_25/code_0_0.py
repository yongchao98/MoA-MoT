def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            decrypted_text.append(chr((ord(char) - shift_amount - shift) % 26 + shift_amount))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "DQOWAZ GB FTQ ZGYNQD AR EOUQZFUEFE MZP MDFUEFE FTMF TQ AR UZHQZFADE FTMF TQ TME GZPQD TUE BMK."

for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")