def decrypt_caesar_cipher(ciphertext, shift):
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

ciphertext = "ESP NZAJCTRSE WLHD ZQ ESP AWLNP HSPCP JZF LCP WZNLEPO LWDZ RZGPCY HSLE JZF NLY OZ HTES ESTD HZCV."

for shift in range(1, 26):
    decrypted_text = decrypt_caesar_cipher(ciphertext, shift)
    print(f"Shift {shift}: {decrypted_text}")