from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_char = [char for char, cnt in count.items() if cnt % 2 != 0]
    
    if len(odd_char) > 1:
        return "Not possible to form a palindrome"
    
    half_palindrome = []
    middle_char = odd_char[0] if odd_char else ''
    
    for char, cnt in count.items():
        half_palindrome.append(char * (cnt // 2))
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + middle_char + half_palindrome[::-1]

letters = ['c', 'f', 'c']
palindrome = form_palindrome(letters)
print(palindrome)