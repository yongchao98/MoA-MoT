from collections import Counter

letters = ['n', 'j', 'n', 'j', 'd', 'j', 's', 's', 'd']
letter_count = Counter(letters)

# Find the center character (odd count)
center = [char for char, count in letter_count.items() if count % 2 == 1][0]

# Construct the half palindrome
half_palindrome = []
for char, count in letter_count.items():
    half_palindrome.extend([char] * (count // 2))

# Form the full palindrome
palindrome = half_palindrome + [center] + half_palindrome[::-1]
print(''.join(palindrome))