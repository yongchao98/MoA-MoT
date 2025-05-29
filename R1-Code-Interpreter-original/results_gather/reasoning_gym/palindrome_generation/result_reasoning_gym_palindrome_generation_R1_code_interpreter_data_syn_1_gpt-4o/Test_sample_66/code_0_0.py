from collections import Counter

# List of letters
letters = ['z', 'k', 'b', 'o', 'o', 'b', 'z', 'k']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_letter = ''

# Determine the half palindrome and middle letter
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_letter = letter
    half_palindrome.extend([letter] * (count // 2))

# Form the palindrome
half_palindrome.sort()
palindrome = ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))

print(palindrome)