from collections import Counter

# Given letters
letters = ['h', 'l', 'f', 'n', 'l', 'n', 'h', 'o', 'f']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_letter = ''

# Construct half of the palindrome and find the middle letter if any
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_letter = letter
    half_palindrome.extend([letter] * (count // 2))

# Sort half_palindrome to ensure consistent output
half_palindrome.sort()

# Construct the full palindrome
full_palindrome = ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))

print(full_palindrome)