from collections import Counter

# Given letters
letters = ['g', 'a', 'g', 'f', 'a', 'f']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Initialize the first half and second half of the palindrome
first_half = []
second_half = []

# Construct the palindrome
for letter, count in letter_count.items():
    half_count = count // 2
    first_half.append(letter * half_count)
    second_half.append(letter * half_count)

# Join the halves to form the palindrome
palindrome = ''.join(first_half) + ''.join(second_half[::-1])

print(palindrome)