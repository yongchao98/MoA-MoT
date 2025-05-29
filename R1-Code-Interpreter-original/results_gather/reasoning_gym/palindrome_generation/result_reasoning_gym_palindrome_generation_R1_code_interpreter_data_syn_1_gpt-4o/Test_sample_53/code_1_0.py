# Define the list of letters
letters = ['v', 's', 'i', 'v', 's', 'i']

# Check if a given string is a palindrome
def is_palindrome(s):
    return s == s[::-1]

# Form a palindrome
palindrome = "visisiv"

# Verify if the formed string is a palindrome
print(is_palindrome(palindrome))