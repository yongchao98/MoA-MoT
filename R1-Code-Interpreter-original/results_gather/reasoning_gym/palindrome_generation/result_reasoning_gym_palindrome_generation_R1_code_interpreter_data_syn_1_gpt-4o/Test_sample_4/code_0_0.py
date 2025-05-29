from collections import Counter

# Given letters
letters = ['x', 'b', 'a', 'b', 'a', 'u', 'x', 'u']

# Count occurrences of each letter
letter_count = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
right_part = []

# Construct the palindrome
for letter, count in letter_count.items():
    # Add half of the count to the left part and half to the right part
    left_part.extend([letter] * (count // 2))
    right_part.extend([letter] * (count // 2))

# The right part should be the reverse of the left part
right_part.reverse()

# Combine the parts to form the palindrome
palindrome = ''.join(left_part + right_part)

print(palindrome)