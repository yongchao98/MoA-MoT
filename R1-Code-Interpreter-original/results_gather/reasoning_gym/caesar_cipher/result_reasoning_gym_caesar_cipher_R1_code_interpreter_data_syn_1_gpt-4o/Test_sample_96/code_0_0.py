def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift % 26
            shifted = ord(char) - shift_amount
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

ciphertext = "FVE QB ABG PBHAG BA ZR SBE URYC"
shift = 13
decrypted_message = caesar_decrypt(ciphertext, shift)
print(decrypted_message)