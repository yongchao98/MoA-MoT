def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
            elif char.islower():
                if shifted < ord('a'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "CWSDR MYXDSXEOC RSC BYEXN KXN OXDOBC DRO RKVV"
for shift in range(1, 26):
    print(f"Shift {shift}: {caesar_decrypt(ciphertext, shift)}")