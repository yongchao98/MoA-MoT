def is_palindrome(s):
    return s == s[::-1]

letters = ['f', 'q', 'y', 'q', 'y', 'f']
palindrome = ''.join(letters)
print(is_palindrome(palindrome))