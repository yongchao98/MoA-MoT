from collections import Counter

# Given letters
letters = ['x', 'x', 's', 'b', 'b', 's', 'm']

# Count occurrences of each letter
letter_count = Counter(letters)

# Separate letters with even and odd counts
even_letters = []
odd_letter = None

for letter, count in letter_count.items():
    if count % 2 == 0:
        even_letters.extend([letter] * (count // 2))
    else:
        odd_letter = letter

# Form the palindrome
half_palindrome = ''.join(even_letters)
if odd_letter:
    palindrome = half_palindrome + odd_letter + half_palindrome[::-1]
else:
    palindrome = half_palindrome + half_palindrome[::-1]

print(palindrome)