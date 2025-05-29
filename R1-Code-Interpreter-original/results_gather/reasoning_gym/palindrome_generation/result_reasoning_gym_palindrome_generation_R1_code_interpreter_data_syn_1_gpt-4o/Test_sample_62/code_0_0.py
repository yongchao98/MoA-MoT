# Define the list of letters
letters = ['d', 'h', 'h', 'n', 'd']

# Attempt to form a palindrome
# We know 'n' can be in the center, and 'd' and 'h' can be mirrored around it
palindrome = "dhnhd"

# Check if the formed string is a palindrome
is_palindrome = palindrome == palindrome[::-1]

# Print the result
print(is_palindrome)