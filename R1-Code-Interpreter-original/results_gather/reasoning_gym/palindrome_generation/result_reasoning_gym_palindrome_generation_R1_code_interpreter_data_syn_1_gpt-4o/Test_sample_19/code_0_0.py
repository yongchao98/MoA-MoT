from collections import Counter

def can_form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count.values() if c % 2 != 0)
    return odd_count <= 1

def form_palindrome(letters):
    if not can_form_palindrome(letters):
        return None
    
    half = []
    middle = []
    
    for letter, cnt in Counter(letters).items():
        half.extend([letter] * (cnt // 2))
        if cnt % 2 == 1:
            middle.append(letter)
    
    return ''.join(half + middle + half[::-1])

letters = ['o', 'o', 'g', 'v', 'g', 'v']
palindrome = form_palindrome(letters)
print(palindrome)