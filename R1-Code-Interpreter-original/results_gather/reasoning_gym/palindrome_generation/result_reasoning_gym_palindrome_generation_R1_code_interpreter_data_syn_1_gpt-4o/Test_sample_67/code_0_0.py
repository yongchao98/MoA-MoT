from collections import Counter

# Given letters
letters = ['c', 'f', 'j', 'd', 'p', 'j', 'f', 'p', 'd']

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

# Form the full palindrome
half_palindrome.sort()
palindrome = ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))

print(palindrome)