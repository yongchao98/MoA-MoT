def is_palindrome(s):
    return s == s[::-1]

letters = ['v', 'v', 'i', 'q', 'q']
palindrome = "qvivq"

print(is_palindrome(palindrome))