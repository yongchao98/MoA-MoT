from collections import Counter

letters = ['h', 't', 't', 't', 'w', 'h', 't', 'w']
letter_count = Counter(letters)

# Construct the first half of the palindrome
half_palindrome = []
for letter, count in letter_count.items():
    half_palindrome.extend([letter] * (count // 2))

# Create the full palindrome
palindrome = half_palindrome + half_palindrome[::-1]

# If there's an odd count, place one in the middle
for letter, count in letter_count.items():
    if count % 2 == 1:
        palindrome.insert(len(half_palindrome), letter)
        break

print(''.join(palindrome))