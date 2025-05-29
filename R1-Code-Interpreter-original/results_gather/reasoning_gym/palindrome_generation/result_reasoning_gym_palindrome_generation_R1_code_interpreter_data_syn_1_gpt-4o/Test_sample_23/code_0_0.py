def is_palindrome(s):
    return s == s[::-1]

letters = ['f', 'z', 'z', 'f']
palindrome = ''.join(letters)
print(is_palindrome(palindrome), palindrome)