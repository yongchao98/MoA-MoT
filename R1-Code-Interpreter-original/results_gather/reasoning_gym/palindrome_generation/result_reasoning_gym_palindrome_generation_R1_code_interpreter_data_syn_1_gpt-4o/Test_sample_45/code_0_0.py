from collections import Counter

def form_palindrome(letters):
    count = Counter(letters)
    odd_char = [char for char, cnt in count.items() if cnt % 2 == 1]
    even_chars = [char for char, cnt in count.items() if cnt % 2 == 0]
    
    # Since we know the input is valid, we can directly form the palindrome
    half_palindrome = []
    for char in even_chars:
        half_palindrome.append(char * (count[char] // 2))
    
    # Add the odd character in the middle if it exists
    middle = odd_char[0] if odd_char else ''
    
    # Form the full palindrome
    half_palindrome = ''.join(half_palindrome)
    full_palindrome = half_palindrome + middle + half_palindrome[::-1]
    
    return full_palindrome

letters = ['v', 'd', 'v', 'g', 'i', 'i', 'g']
palindrome = form_palindrome(letters)
print(palindrome)