def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_char = chr((ord(char) - start - shift_amount) % 26 + start)
            decrypted_text.append(decrypted_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

encrypted_text = "WNA KLAJEJC PDA DA ATLHWEJAZ"
for shift in range(1, 26):
    decrypted_text = decrypt_caesar_cipher(encrypted_text, shift)
    print(f"Shift {shift}: {decrypted_text}")