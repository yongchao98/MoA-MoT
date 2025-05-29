def is_palindrome(s):
    return s == s[::-1]

letters = ['a', 'a', 'c', 'c', 'p', 'p', 'n', 'n']
palindrome = "acpnnpca"

print(is_palindrome(palindrome))