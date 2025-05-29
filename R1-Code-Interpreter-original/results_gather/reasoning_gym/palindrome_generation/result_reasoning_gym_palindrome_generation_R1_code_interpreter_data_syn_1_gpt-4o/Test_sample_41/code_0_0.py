from collections import Counter

# Given letters
letters = ['e', 'd', 'j', 'd', 'e', 'm', 'm']

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
half_palindrome.sort()
palindrome = ''.join(half_palindrome) + center + ''.join(reversed(half_palindrome))

print(palindrome)