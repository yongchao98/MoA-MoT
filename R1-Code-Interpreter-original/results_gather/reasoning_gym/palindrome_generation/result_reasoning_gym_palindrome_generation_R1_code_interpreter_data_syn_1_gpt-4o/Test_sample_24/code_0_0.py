from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_char = [char for char, cnt in count.items() if cnt % 2 == 1]
    even_chars = [char for char, cnt in count.items() if cnt % 2 == 0]
    
    # Start forming the palindrome
    half_palindrome = []
    for char, cnt in count.items():
        half_palindrome.append(char * (cnt // 2))
    
    # Join the half palindrome
    half_palindrome = ''.join(half_palindrome)
    
    # If there's an odd character, it goes in the middle
    if odd_char:
        palindrome = half_palindrome + odd_char[0] + half_palindrome[::-1]
    else:
        palindrome = half_palindrome + half_palindrome[::-1]
    
    return palindrome

letters = ['q', 'm', 'm']
palindrome = form_palindrome(letters)
print(palindrome)