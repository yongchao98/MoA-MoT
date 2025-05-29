from collections import Counter

# Given letters
letters = ['l', 'b', 'l', 's', 's', 'b', 's', 's', 'l']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
center = ''

# Construct half of the palindrome and find the center
for letter, count in letter_count.items():
    if count % 2 == 1:
        center = letter
    half_palindrome.extend([letter] * (count // 2))

# Form the full palindrome
full_palindrome = ''.join(half_palindrome) + center + ''.join(reversed(half_palindrome))

print(full_palindrome)