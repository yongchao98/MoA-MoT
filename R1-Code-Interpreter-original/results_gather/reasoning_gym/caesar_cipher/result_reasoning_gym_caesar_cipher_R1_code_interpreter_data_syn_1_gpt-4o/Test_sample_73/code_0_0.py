def decrypt_caesar_cipher(text, shift):
    decrypted_text = []
    for char in text:
        if char.isalpha():
            shift_amount = shift if char.isupper() else shift
            new_char = chr((ord(char) - shift_amount - 65) % 26 + 65)
            decrypted_text.append(new_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

cipher_text = "FTQ QUZ AD RQPQDMX FMJ UPQZFURUOMFUAZ ZGYNQD UE"

for shift in range(1, 26):
    decrypted = decrypt_caesar_cipher(cipher_text, shift)
    print(f"Shift {shift}: {decrypted}")