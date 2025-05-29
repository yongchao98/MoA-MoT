from collections import Counter

def form_palindrome(letters):
    letter_count = Counter(letters)
    half_palindrome = []
    
    for letter, count in letter_count.items():
        half_palindrome.append(letter * (count // 2))
    
    half_palindrome = ''.join(half_palindrome)
    palindrome = half_palindrome + half_palindrome[::-1]
    
    return palindrome

letters = ['h', 'a', 'h', 'a']
palindrome = form_palindrome(letters)
print(palindrome)