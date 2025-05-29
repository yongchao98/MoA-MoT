from collections import Counter

# List of letters
letters = ['c', 's', 'w', 'd', 'x', 'x', 'd', 'w', 's', 'c']

# Count the occurrences of each letter
letter_counts = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
right_part = []

# Construct the palindrome
for letter, count in letter_counts.items():
    # Add half of the count to the left part and half to the right part
    left_part.extend([letter] * (count // 2))
    right_part.extend([letter] * (count // 2))

# The palindrome is the left part + right part reversed
palindrome = left_part + right_part[::-1]

# Convert the list to a string
palindrome_str = ''.join(palindrome)

# Print the result
print(palindrome_str)