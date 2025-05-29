def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift if char.isupper() else shift
            new_char = chr((ord(char) - shift_amount - 65) % 26 + 65) if char.isupper() else chr((ord(char) - shift_amount - 97) % 26 + 97)
            decrypted_text.append(new_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "X AFPZRPPFLK TXP DLFKD LK TEBK EB BKQBOBA"
for shift in range(1, 26):
    print(f"Shift {shift}: {caesar_decrypt(ciphertext, shift)}")