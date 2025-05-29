def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shifted = ord(char) - shift
            if char.isupper():
                if shifted < ord('A'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
            else:
                if shifted < ord('a'):
                    shifted += 26
                decrypted_text.append(chr(shifted))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "QRB RAXW LXWBCRCDCRXW RB CJGNM CX CQN DCVXBC KH CQN QNJEH BCAJRW CQJC RB YDC DYXW RC."

for shift in range(1, 26):
    decrypted_text = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_text}")