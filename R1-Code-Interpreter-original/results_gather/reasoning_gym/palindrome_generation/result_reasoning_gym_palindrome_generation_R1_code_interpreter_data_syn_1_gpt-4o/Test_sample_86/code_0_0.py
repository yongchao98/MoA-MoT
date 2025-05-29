from collections import Counter

# Given letters
letters = ['i', 'j', 's', 'i', 's']

# Count the frequency of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_char = ''

# Determine the half palindrome and middle character
for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.extend([char] * (count // 2))

# Form the palindrome
first_half = ''.join(half_palindrome)
second_half = first_half[::-1]
palindrome = first_half + middle_char + second_half

print(palindrome)