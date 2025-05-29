letters = ['t', 'r', 'w', 'w', 't']

# Count occurrences of each letter
from collections import Counter

letter_count = Counter(letters)

# Construct the palindrome
half_palindrome = []
middle_char = ''

# Find the middle character (if any) and construct half of the palindrome
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_char = letter
    half_palindrome.extend([letter] * (count // 2))

# Form the full palindrome
palindrome = ''.join(half_palindrome) + middle_char + ''.join(reversed(half_palindrome))

print(palindrome)