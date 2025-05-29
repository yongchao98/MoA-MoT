from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    middle_letter = None

    for letter, count in letter_count.items():
        if count % 2 != 0:
            if middle_letter is None:
                middle_letter = letter
            else:
                return None  # More than one letter with an odd count, can't form a palindrome
        half_palindrome.extend([letter] * (count // 2))

    # Form the full palindrome
    full_palindrome = half_palindrome + ([middle_letter] if middle_letter else []) + half_palindrome[::-1]
    return ''.join(full_palindrome)

letters = ['y', 'h', 'd', 's', 'y', 'h', 'd', 's']
palindrome = form_palindrome(letters)
print(palindrome)