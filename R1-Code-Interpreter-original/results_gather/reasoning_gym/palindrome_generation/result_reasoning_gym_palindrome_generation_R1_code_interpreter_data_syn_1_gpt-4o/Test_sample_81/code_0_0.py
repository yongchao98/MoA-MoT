from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count if count[c] % 2 != 0)
    
    if odd_count > 1:
        return "Not possible to form a palindrome"
    
    half_palindrome = []
    middle = []
    
    for char, cnt in count.items():
        if cnt % 2 != 0:
            middle.append(char)
        half_palindrome.append(char * (cnt // 2))
    
    half_palindrome = ''.join(half_palindrome)
    middle = ''.join(middle)
    
    return half_palindrome + middle + half_palindrome[::-1]

letters = ['r', 'z', 'r']
palindrome = form_palindrome(letters)
print(palindrome)