from collections import Counter

# Given letters
letters = ['f', 'c', 'c', 'c', 'y', 'c', 'f']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_letter = ''

# Construct the half palindrome and find the middle letter if any
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_letter = letter
    half_palindrome.append(letter * (count // 2))

# Form the full palindrome
half_palindrome = ''.join(half_palindrome)
palindrome = half_palindrome + middle_letter + half_palindrome[::-1]

print(palindrome)