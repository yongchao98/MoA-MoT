def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
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

ciphertext = "KLV VSHFLDO SURMHFW KDG WR GR ZLWK WKH FLWB RI DV HYHUBERGB VRPH ILIWHHQ PLOHV LQODQG."

for shift in range(1, 26):
    decrypted = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted}")