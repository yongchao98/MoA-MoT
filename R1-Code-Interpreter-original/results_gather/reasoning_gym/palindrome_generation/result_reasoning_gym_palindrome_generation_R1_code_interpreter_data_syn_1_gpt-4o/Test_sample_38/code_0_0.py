def is_palindrome(s):
    return s == s[::-1]

letters = ['l', 'h', 'h', 'l']
# Try to form a palindrome
palindrome = ''.join(letters)
if is_palindrome(palindrome):
    print(palindrome)
else:
    # Try another arrangement
    palindrome = ''.join(letters[::-1])
    if is_palindrome(palindrome):
        print(palindrome)