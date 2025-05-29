from collections import Counter

def can_form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count if count[c] % 2 != 0)
    return odd_count <= 1

def form_palindrome(letters):
    if not can_form_palindrome(letters):
        return None
    
    count = Counter(letters)
    half_palindrome = []
    middle = []
    
    for char, freq in count.items():
        half_palindrome.append(char * (freq // 2))
        if freq % 2 == 1:
            middle.append(char)
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + ''.join(middle) + half_palindrome[::-1]

letters = ['p', 'y', 'h', 'y', 'p', 'h']
palindrome = form_palindrome(letters)
print(palindrome)