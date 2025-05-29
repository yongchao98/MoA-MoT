from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count if count[c] % 2 != 0)
    
    if odd_count > 1:
        return None
    
    half_palindrome = []
    middle = ''
    
    for char, cnt in count.items():
        if cnt % 2 != 0:
            middle = char
        half_palindrome.append(char * (cnt // 2))
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + middle + half_palindrome[::-1]

letters = ['o', 'o', 'k']
palindrome = form_palindrome(letters)
print(palindrome)