def is_palindrome(s):
    return s == s[::-1]

letters = ['f', 'y', 'q', 'q', 'y', 'f']
palindrome = ''.join(letters)
print(is_palindrome(palindrome))