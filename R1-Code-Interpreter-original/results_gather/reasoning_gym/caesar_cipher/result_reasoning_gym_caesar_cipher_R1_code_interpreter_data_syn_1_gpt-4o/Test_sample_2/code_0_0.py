def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_text.append(chr((ord(char) - start - shift_amount) % 26 + start))
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "MRFQD M RQI EQOAZPE M RMUZF DGYNXUZS IME ITUOT SDQI XAGPQD MZP XAGPQD"

for shift in range(1, 26):
    decrypted = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted}")