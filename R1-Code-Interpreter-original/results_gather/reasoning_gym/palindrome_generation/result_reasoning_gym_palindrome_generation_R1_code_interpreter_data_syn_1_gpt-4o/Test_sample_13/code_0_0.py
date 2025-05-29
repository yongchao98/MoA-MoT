from collections import Counter

# Given letters
letters = ['k', 'a', 'k', 'h', 'g', 'g', 'h']

# Count the frequency of each letter
letter_count = Counter(letters)

# Separate letters with even and odd frequencies
even_part = []
odd_part = []

for letter, count in letter_count.items():
    if count % 2 == 0:
        even_part.extend([letter] * (count // 2))
    else:
        odd_part.append(letter)

# Form the palindrome
palindrome = even_part + odd_part + even_part[::-1]

# Join the list to form a string
palindrome_str = ''.join(palindrome)

print(palindrome_str)