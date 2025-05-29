from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    middle_letter = None

    for letter, count in letter_count.items():
        if count % 2 != 0:
            middle_letter = letter
        half_palindrome.extend([letter] * (count // 2))

    # Form the full palindrome
    full_palindrome = half_palindrome + ([middle_letter] if middle_letter else []) + half_palindrome[::-1]
    return ''.join(full_palindrome)

letters = ['n', 'k', 'a', 'k', 'j', 'n', 'a', 'j']
palindrome = form_palindrome(letters)
print(palindrome)