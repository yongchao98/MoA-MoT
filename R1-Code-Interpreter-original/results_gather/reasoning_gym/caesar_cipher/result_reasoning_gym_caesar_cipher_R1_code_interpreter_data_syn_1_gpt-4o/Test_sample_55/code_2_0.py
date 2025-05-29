def calculate_shift(most_frequent_cipher_letter, most_frequent_plain_letter='E'):
    shift = (ord(most_frequent_cipher_letter) - ord(most_frequent_plain_letter)) % 26
    return shift

most_frequent_cipher_letter = 'C'
shift = calculate_shift(most_frequent_cipher_letter)
print(shift)