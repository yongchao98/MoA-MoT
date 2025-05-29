from collections import Counter

# Given letters
letters = ['r', 'g', 'q', 'n', 'g', 'p', 'p', 'r', 'q', 'n']

# Count occurrences of each letter
letter_count = Counter(letters)

# Separate letters into even and odd counts
odd_count_letters = [letter for letter, count in letter_count.items() if count % 2 != 0]
even_count_letters = [letter for letter, count in letter_count.items() if count % 2 == 0]

# Check if a palindrome can be formed
if len(odd_count_letters) > 1:
    print("Cannot form a palindrome")
else:
    # Construct the palindrome
    half_palindrome = []
    for letter, count in letter_count.items():
        half_palindrome.append(letter * (count // 2))
    
    # Sort to ensure consistent output
    half_palindrome.sort()
    
    # Form the full palindrome
    if odd_count_letters:
        center_letter = odd_count_letters[0]
        palindrome = ''.join(half_palindrome) + center_letter + ''.join(half_palindrome[::-1])
    else:
        palindrome = ''.join(half_palindrome) + ''.join(half_palindrome[::-1])
    
    print(palindrome)