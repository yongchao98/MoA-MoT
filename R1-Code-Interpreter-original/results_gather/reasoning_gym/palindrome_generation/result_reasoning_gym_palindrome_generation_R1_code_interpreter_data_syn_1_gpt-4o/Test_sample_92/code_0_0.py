from collections import Counter

letters = ['z', 'o', 'o', 'e', 'z', 'e', 'n']
letter_count = Counter(letters)

# Find the middle character (if any)
middle_char = ''
half_palindrome = []

for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.extend([char] * (count // 2))

# Sort half_palindrome to ensure consistent output
half_palindrome.sort()

# Form the palindrome
palindrome = ''.join(half_palindrome) + middle_char + ''.join(reversed(half_palindrome))
print(palindrome)