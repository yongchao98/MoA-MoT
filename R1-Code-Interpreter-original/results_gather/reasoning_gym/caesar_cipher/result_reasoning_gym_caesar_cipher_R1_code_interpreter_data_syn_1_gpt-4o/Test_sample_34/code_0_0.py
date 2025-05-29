def caesar_decrypt(ciphertext, shift):
    decrypted_text = []
    for char in ciphertext:
        if char.isalpha():
            shift_amount = shift
            if char.isupper():
                start = ord('A')
            else:
                start = ord('a')
            decrypted_char = chr(start + (ord(char) - start - shift_amount) % 26)
            decrypted_text.append(decrypted_char)
        else:
            decrypted_text.append(char)
    return ''.join(decrypted_text)

ciphertext = "CWSDR NYOC SC DY MYXXOMD RSC DRO GSBOC YP GRSMR MYWWEXSMKDO GSDR RSC ZKBSC WKXCSYX."

# Try all possible shifts
for shift in range(1, 26):
    decrypted_text = caesar_decrypt(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_text}")