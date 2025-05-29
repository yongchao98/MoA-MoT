from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    half_palindrome = []
    
    for char, freq in count.items():
        half_palindrome.append(char * (freq // 2))
    
    half_palindrome = ''.join(half_palindrome)
    return half_palindrome + half_palindrome[::-1]

letters = ['s', 't', 's', 't']
palindrome = form_palindrome(letters)
print(palindrome)