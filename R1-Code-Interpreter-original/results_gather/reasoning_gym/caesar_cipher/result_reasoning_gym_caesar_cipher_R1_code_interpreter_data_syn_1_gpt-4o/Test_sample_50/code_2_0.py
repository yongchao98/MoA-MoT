def rot13(text):
    result = []
    for char in text:
        if char.isalpha():
            shift_amount = 65 if char.isupper() else 97
            result.append(chr((ord(char) - shift_amount + 13) % 26 + shift_amount))
        else:
            result.append(char)
    return ''.join(result)

ciphertext = "IUFT FTQUD TQMPE NQZF PAIZ AHQD FTQUD QXQOFDUO FTUDFK EOUQZFURUO YQZ IQDQ MNEADOQP UZ FDMZEOQZPQZFMX OMXOGXMFUAZE"
decrypted = rot13(ciphertext)
print(decrypted)