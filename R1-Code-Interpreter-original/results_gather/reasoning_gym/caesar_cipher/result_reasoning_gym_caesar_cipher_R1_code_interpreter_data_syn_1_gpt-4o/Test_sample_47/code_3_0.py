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

rot13_text = "JVGU GURVE URNQF ORAG QBJA BIRE GURVE RYRPGEVP GUVEGL FPVRAGVSVP ZRA JRER NOFBEORQ VA GENAFPRAQRAGNY PNYPHYNGVBAF."

# Decrypting the ROT13 text
decrypted = caesar_decrypt(rot13_text, 13)
print(decrypted)