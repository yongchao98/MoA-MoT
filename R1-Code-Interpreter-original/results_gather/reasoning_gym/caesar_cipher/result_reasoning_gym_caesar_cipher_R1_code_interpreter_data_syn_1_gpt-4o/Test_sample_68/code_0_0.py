def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            decrypted_text.append(chr((ord(char) - shift_amount - shift) % 26 + shift_amount))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "LXIW IWTXG WTPSH QTCI SDLC DKTG IWTXG TATRIGXR IWXGIN HRXTCIXUXR BTC LTGT PQHDGQTS XC IGPCHRTCSTCIPA RPARJAPIXDCH."

for shift in range(1, 26):
    decrypted = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted}")