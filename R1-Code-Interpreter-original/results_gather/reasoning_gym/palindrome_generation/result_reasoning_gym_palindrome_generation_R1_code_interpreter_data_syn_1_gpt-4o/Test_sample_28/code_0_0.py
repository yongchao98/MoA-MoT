from collections import Counter

letters = ['p', 'f', 'h', 'f', 'e', 'h', 'e']
letter_count = Counter(letters)

# Find the middle character (if any)
middle_char = [char for char, count in letter_count.items() if count % 2 == 1]

# Construct the first half of the palindrome
half_palindrome = []
for char, count in letter_count.items():
    half_palindrome.append(char * (count // 2))

# Join the first half, middle character, and the reverse of the first half
first_half = ''.join(half_palindrome)
palindrome = first_half + ''.join(middle_char) + first_half[::-1]

print(palindrome)