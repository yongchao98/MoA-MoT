def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            decrypted_text.append(chr((ord(char) - shift_amount - shift) % 26 + shift_amount))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "UFERKZFEJ RIV RTTVGKVU ZE R ELDSVI FW FKYVI NRPJ ZETCLUZEX FECZEV GRPDVEKJ REU TIVUZK TRIU UFERKZFEJ."

for shift in range(10, 26):
    decrypted_text = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_text}")