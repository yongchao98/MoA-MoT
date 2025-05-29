from collections import Counter

def form_palindrome(letters):
    # Count the frequency of each letter
    letter_count = Counter(letters)
    
    # Create two halves of the palindrome
    half_palindrome = []
    
    # Add half of each letter's count to the half_palindrome
    for letter, count in letter_count.items():
        half_palindrome.append(letter * (count // 2))
    
    # Join the half_palindrome to form the first half
    first_half = ''.join(half_palindrome)
    
    # The second half is the reverse of the first half
    second_half = first_half[::-1]
    
    # Form the full palindrome
    palindrome = first_half + second_half
    
    return palindrome

letters = ['f', 'v', 'h', 'f', 'v', 'h']
palindrome = form_palindrome(letters)
print(palindrome)