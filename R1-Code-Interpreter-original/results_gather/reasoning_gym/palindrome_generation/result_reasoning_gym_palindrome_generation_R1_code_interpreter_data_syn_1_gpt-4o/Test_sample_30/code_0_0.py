# Given letters
letters = ['o', 'b', 'b', 'v', 'v']

# Count the frequency of each letter
from collections import Counter

letter_count = Counter(letters)

# Construct the palindrome
half_palindrome = []
middle_char = ''

# Find the middle character and half of the palindrome
for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.extend([char] * (count // 2))

# Form the full palindrome
palindrome = ''.join(half_palindrome) + middle_char + ''.join(reversed(half_palindrome))

print(palindrome)