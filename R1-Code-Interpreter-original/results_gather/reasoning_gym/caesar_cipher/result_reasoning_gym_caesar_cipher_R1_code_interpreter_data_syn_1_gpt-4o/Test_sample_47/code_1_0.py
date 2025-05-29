def rot13(text):
    result = []
    for char in text:
        if char.isalpha():
            shift = 13
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            result.append(chr((ord(char) - start + shift) % 26 + start))
        else:
            result.append(char)
    return ''.join(result)

encoded_text = "JVGU GURVE URNQF ORAG QBJA BIRE GURVE RYRPGEVP GUVEGL FPVRAGVSVP ZRA JRER NOFBEORQ VA GENAFPRAQRAGNY PNYPHYNGVBAF."
decrypted_text = rot13(encoded_text)
print(decrypted_text)