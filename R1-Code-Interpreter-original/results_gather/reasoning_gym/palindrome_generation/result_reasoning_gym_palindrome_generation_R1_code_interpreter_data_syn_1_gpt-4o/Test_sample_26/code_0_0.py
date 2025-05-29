from collections import Counter

def can_form_palindrome(letters):
    letter_count = Counter(letters)
    odd_count = sum(1 for count in letter_count.values() if count % 2 != 0)
    return odd_count <= 1

def form_palindrome(letters):
    if not can_form_palindrome(letters):
        return None
    
    letter_count = Counter(letters)
    half_palindrome = []
    middle_letter = ''
    
    for letter, count in letter_count.items():
        if count % 2 != 0:
            middle_letter = letter
        half_palindrome.extend([letter] * (count // 2))
    
    return ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))

letters = ['f', 'q', 'y', 'q', 'y', 'f']
palindrome = form_palindrome(letters)
print(palindrome)