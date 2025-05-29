from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    
    for letter, count in letter_count.items():
        half_palindrome.append(letter * (count // 2))
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + half_palindrome[::-1]

letters = ['y', 'i', 'i', 'y']
palindrome = form_palindrome(letters)
print(palindrome)