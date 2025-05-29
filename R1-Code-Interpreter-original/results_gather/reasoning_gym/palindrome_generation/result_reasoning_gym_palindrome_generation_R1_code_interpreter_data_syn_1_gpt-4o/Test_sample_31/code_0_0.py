from collections import Counter

# List of letters
letters = ['y', 'l', 'r', 'g', 'l', 'r', 'k', 'g', 'k', 'y']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
right_part = []

# Construct the palindrome
for letter, count in letter_count.items():
    # Add half of the count to the left and half to the right
    left_part.append(letter * (count // 2))
    right_part.append(letter * (count // 2))

# Join the parts to form the palindrome
palindrome = ''.join(left_part) + ''.join(right_part[::-1])

print(palindrome)