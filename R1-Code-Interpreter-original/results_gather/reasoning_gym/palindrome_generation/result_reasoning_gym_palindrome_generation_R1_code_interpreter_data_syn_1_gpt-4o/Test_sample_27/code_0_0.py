from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    
    for letter, count in letter_count.items():
        half_palindrome.extend([letter] * (count // 2))
    
    # Form the full palindrome
    full_palindrome = half_palindrome + half_palindrome[::-1]
    return ''.join(full_palindrome)

letters = ['q', 'l', 'c', 'q', 'l', 'c']
palindrome = form_palindrome(letters)
print(palindrome)