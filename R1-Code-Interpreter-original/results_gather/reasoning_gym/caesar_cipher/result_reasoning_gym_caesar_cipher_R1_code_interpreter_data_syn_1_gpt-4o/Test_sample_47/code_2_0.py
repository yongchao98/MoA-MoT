def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_text.append(chr((ord(char) - start - shift_amount) % 26 + start))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "LXIW IWTXG WTPSH QTCI SDLC DKTG IWTXG TATRIGXR IWXGIN HRXTCIXUXR BTC LTGT PQHDGQTS XC IGPCHRTCSTCIPA RPARJAPIXDCH."

# Trying shift 13 to see if it results in a meaningful sentence
decrypted = caesar_decrypt(ciphertext, 13)
print(f"Shift 13: {decrypted}")