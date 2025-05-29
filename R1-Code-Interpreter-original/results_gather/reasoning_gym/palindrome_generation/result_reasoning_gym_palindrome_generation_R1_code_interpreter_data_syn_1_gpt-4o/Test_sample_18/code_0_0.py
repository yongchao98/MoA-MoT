from collections import Counter

# Given letters
letters = ['t', 'j', 'w', 'j', 'x', 't', 'w', 's', 's', 'x']

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

# The palindrome is the left part + middle (if any) + reversed right part
palindrome = left_part + right_part[::-1]

# Join the list into a string
palindrome_str = ''.join(palindrome)

print(palindrome_str)