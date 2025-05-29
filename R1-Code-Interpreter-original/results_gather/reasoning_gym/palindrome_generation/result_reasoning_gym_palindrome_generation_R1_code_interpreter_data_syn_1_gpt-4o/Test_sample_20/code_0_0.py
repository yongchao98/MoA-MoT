from collections import Counter

letters = ['l', 'b', 's', 'm', 'b', 's', 'm']
letter_count = Counter(letters)

# Find the middle character (if any)
middle_char = [char for char, count in letter_count.items() if count % 2 == 1]

# Prepare the half part of the palindrome
half_palindrome = []
for char, count in letter_count.items():
    half_palindrome.extend([char] * (count // 2))

# Construct the full palindrome
half_palindrome.sort()
palindrome = half_palindrome + middle_char + half_palindrome[::-1]

# Join to form the final palindrome string
palindrome_str = ''.join(palindrome)
print(palindrome_str)