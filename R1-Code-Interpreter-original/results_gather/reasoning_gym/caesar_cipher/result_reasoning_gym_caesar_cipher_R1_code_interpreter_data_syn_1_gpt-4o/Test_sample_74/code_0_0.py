def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
            elif char.islower():
                if shifted < ord('a'):
                    shifted += 26
            decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "UIF DBTLFU CFJOH MBJE EPXO JO UIF NJEEMF PG UIF UIF UFMFQIPUF XBT HPU JO SFBEJOFTT."
shift = 1
decrypted_text = decrypt_caesar_cipher(cipher_text, shift)
print(decrypted_text)