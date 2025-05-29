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

cipher_text = "ESP NZAJCTRSE WLHD ZQ ESP AWLNP HSPCP JZF LCP WZNLEPO LWDZ RZGPCY HSLE JZF NLY OZ HTES ESTD HZCV."
shift = 18
decrypted_message = decrypt_caesar_cipher(cipher_text, shift)
print(decrypted_message)