from collections import Counter

# Given letters
letters = ['r', 'j', 'j', 'r', 'r']

# Count the frequency of each letter
letter_count = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
middle_part = []

# Construct the palindrome
for letter, count in letter_count.items():
    if count % 2 == 1:
        middle_part.append(letter * (count % 2))
    left_part.append(letter * (count // 2))

# Form the palindrome
palindrome = ''.join(left_part) + ''.join(middle_part) + ''.join(reversed(left_part))
print(palindrome)