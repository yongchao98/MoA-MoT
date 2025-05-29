def decrypt_caesar_cipher(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift if char.isupper() else shift
            new_char = chr((ord(char) - shift_amount - 65) % 26 + 65)
            decrypted_text.append(new_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "FUR VF URE URNQ FHAX VA GUR CVYYBJF."

for shift in range(1, 26):
    decrypted_message = decrypt_caesar_cipher(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_message}")