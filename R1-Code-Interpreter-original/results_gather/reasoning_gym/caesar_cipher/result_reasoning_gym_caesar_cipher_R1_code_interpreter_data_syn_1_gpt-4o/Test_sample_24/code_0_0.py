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

cipher_text = "BOE BOZ BEEJUJPOBM UFSNT JNQPTFE CZ UIF DPQZSJHIU IPMEFS."

# Try all possible shifts
for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")