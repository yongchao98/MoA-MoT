from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    half_palindrome = []
    
    for letter, freq in count.items():
        half_palindrome.append(letter * (freq // 2))
    
    # Join the half palindrome and its reverse to form the full palindrome
    half_palindrome = ''.join(half_palindrome)
    full_palindrome = half_palindrome + half_palindrome[::-1]
    
    return full_palindrome

letters = ['s', 'q', 's', 'q', 't', 't']
palindrome = form_palindrome(letters)
print(palindrome)