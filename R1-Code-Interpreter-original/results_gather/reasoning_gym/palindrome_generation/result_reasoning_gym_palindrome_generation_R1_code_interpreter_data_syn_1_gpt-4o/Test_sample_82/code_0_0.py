from collections import Counter

# Given letters
letters = ['o', 'i', 'o', 'z', 'i', 'z']

# Count the frequency of each letter
letter_count = Counter(letters)

# Initialize parts of the palindrome
half_palindrome = []
middle_char = ''

# Determine the half palindrome and middle character
for char, count in letter_count.items():
    if count % 2 == 1:
        middle_char = char
    half_palindrome.append(char * (count // 2))

# Form the first half of the palindrome
first_half = ''.join(half_palindrome)

# Form the full palindrome
palindrome = first_half + middle_char + first_half[::-1]

# Output the result
print(palindrome)