from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    middle_letter = None
    
    for letter, count in letter_count.items():
        if count % 2 == 1:
            middle_letter = letter
        half_palindrome.extend([letter] * (count // 2))
    
    half_palindrome.sort()
    if middle_letter:
        return ''.join(half_palindrome) + middle_letter + ''.join(reversed(half_palindrome))
    else:
        return ''.join(half_palindrome) + ''.join(reversed(half_palindrome))

letters = ['x', 'o', 'o', 'j', 'x', 'j']
palindrome = form_palindrome(letters)
print(palindrome)