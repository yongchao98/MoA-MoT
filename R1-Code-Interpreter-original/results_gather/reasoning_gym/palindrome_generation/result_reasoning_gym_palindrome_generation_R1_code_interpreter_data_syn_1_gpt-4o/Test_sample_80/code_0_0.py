from collections import Counter

# Given letters
letters = ['h', 'k', 'm', 'k', 'h', 'v', 'p', 'v', 'm', 'p']

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

# The right part should be reversed to form a palindrome
right_part.reverse()

# Combine the parts to form the palindrome
palindrome = ''.join(left_part + right_part)

print(palindrome)