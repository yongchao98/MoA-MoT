from collections import Counter

# Given letters
letters = ['i', 'i', 'd', 'm', 'm', 'a', 'n', 'a', 'n', 'd']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Initialize the left and right parts of the palindrome
left_part = []
right_part = []

# Construct the palindrome
for letter, count in letter_count.items():
    half_count = count // 2
    left_part.append(letter * half_count)
    right_part.append(letter * half_count)

# Join the parts to form the palindrome
palindrome = ''.join(left_part) + ''.join(right_part[::-1])

print(palindrome)