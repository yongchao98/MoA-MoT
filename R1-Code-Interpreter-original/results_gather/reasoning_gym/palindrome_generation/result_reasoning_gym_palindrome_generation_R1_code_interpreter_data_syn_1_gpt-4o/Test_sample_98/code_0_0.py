from collections import Counter

# Given letters
letters = ['q', 'q', 'a', 'u', 't', 'h', 't', 'h', 'u']

# Count occurrences of each letter
letter_counts = Counter(letters)

# Initialize parts of the palindrome
left_half = []
center = ""

# Construct the palindrome
for letter, count in letter_counts.items():
    if count % 2 == 1:
        center = letter
    left_half.extend([letter] * (count // 2))

# Form the palindrome
right_half = left_half[::-1]
palindrome = ''.join(left_half) + center + ''.join(right_half)

print(palindrome)