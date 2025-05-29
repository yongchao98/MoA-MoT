def is_palindrome(s):
    return s == s[::-1]

letters = ['d', 'k', 'k', 'd']
palindrome = ''.join(letters)

print(is_palindrome(palindrome), palindrome)