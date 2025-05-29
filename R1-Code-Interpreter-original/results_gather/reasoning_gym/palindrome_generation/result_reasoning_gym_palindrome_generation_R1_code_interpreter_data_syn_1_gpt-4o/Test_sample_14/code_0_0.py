def is_palindrome(s):
    return s == s[::-1]

letters = ['p', 'p', 'k', 'k']
# Try to form a palindrome
palindrome = "pkkp"

# Check if the formed string is a palindrome
print(is_palindrome(palindrome))