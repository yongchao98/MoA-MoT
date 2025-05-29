from collections import Counter

letters = ['b', 'h', 'b', 'f', 'l', 'l', 'f']
letter_count = Counter(letters)

# Find the middle character (if any)
middle_char = [char for char, count in letter_count.items() if count % 2 == 1]

# Find pairs of characters
pairs = [char * (count // 2) for char, count in letter_count.items() if count // 2 > 0]

# Construct the palindrome
half_palindrome = ''.join(pairs)
palindrome = half_palindrome + ''.join(middle_char) + half_palindrome[::-1]

print(palindrome)