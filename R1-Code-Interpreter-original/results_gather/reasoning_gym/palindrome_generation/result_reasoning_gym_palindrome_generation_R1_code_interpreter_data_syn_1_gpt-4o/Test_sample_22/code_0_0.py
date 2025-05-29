from collections import Counter

# Given letters
letters = ['h', 'h', 'l', 'p', 'l']

# Count the frequency of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_char = ''

# Construct half of the palindrome and find the middle character if needed
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_char = letter
    half_palindrome.append(letter * (count // 2))

# Form the full palindrome
half_palindrome = ''.join(half_palindrome)
palindrome = half_palindrome + middle_char + half_palindrome[::-1]

print(palindrome)