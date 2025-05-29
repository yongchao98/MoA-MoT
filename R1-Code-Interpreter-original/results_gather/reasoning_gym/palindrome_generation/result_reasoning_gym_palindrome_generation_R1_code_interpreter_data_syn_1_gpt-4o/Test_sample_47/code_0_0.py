from collections import Counter

# Given letters
letters = ['p', 'p', 'l', 'h', 'l']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_char = ''

# Construct the half palindrome and find the middle character if needed
for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.extend([char] * (count // 2))

# Form the full palindrome
palindrome = ''.join(half_palindrome) + middle_char + ''.join(reversed(half_palindrome))

print(palindrome)