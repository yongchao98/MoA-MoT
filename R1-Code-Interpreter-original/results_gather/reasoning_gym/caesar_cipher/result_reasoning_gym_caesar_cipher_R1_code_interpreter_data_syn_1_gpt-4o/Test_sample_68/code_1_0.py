def rot13_decrypt(text):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            decrypted_text.append(chr((ord(char) - shift_amount + 13) % 26 + shift_amount))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

encrypted_text = "JVGU GURVE URNQF ORAG QBJA BIRE GURVE RYRPGEVP GUVEGL FPVRAGVSVP ZRA JRER NOFBEORQ VA GENAFPRAQRAGNY PNYPHYNGVBAF"
decrypted_text = rot13_decrypt(encrypted_text)
print(decrypted_text)