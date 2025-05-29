from collections import Counter

# Given letters
letters = ['w', 'i', 'k', 'k', 'z', 'i', 'x', 'x', 'z']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_letter = ''

# Construct the half palindrome and find the middle letter if any
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_letter = letter
    half_palindrome.extend([letter] * (count // 2))

# Sort the half palindrome for consistency
half_palindrome.sort()

# Form the full palindrome
full_palindrome = ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))

print(full_palindrome)