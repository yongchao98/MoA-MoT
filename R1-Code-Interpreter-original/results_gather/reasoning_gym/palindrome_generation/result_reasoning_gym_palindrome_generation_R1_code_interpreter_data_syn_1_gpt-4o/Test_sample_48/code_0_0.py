from collections import Counter

# Given letters
letters = ['w', 'k', 'w', 'n', 'x', 'k', 'x']

# Count the occurrences of each letter
letter_count = Counter(letters)

# Separate letters with even and odd counts
even_part = []
odd_part = []

for letter, count in letter_count.items():
    if count % 2 == 0:
        even_part.extend([letter] * (count // 2))
    else:
        odd_part.append(letter)

# Construct the palindrome
# Since we know there's only one odd, it will be in the middle
middle = odd_part[0] if odd_part else ''
half_palindrome = ''.join(even_part)

# Form the full palindrome
palindrome = half_palindrome + middle + half_palindrome[::-1]

print(palindrome)