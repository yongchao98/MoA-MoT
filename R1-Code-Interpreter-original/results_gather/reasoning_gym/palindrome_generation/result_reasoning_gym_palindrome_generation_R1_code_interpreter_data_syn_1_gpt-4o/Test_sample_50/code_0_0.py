from collections import Counter

# Given letters
letters = ['z', 'u', 'j', 'z', 'u']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_char = ''

# Construct half of the palindrome and find the middle character if needed
for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.extend([char] * (count // 2))

# Form the full palindrome
palindrome = ''.join(half_palindrome) + middle_char + ''.join(reversed(half_palindrome))

print(palindrome)