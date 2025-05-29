from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_count = sum(1 for c in count if count[c] % 2 != 0)
    
    if odd_count > 1:
        return "Not possible to form a palindrome"
    
    half_palindrome = []
    middle = ""
    
    for letter, freq in count.items():
        if freq % 2 != 0:
            middle = letter
        half_palindrome.append(letter * (freq // 2))
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + middle + half_palindrome[::-1]

letters = ['y', 'g', 'g']
result = form_palindrome(letters)
print(result)