from collections import Counter

# Given letters
letters = ['c', 'x', 'p', 'p', 'y', 'x', 'c', 't', 't', 'y']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
middle_part = []

# Construct the palindrome
for letter, count in letter_count.items():
    # Add half of the count to the left part
    left_part.extend([letter] * (count // 2))
    # If there's an odd count, it can be a middle part
    if count % 2 == 1:
        middle_part.append(letter)

# The right part is the reverse of the left part
right_part = left_part[::-1]

# Form the palindrome
palindrome = ''.join(left_part + middle_part + right_part)

print(palindrome)