import string

# Known incorrect letters
incorrect_letters = set("FYKSQBIRULEDNMW")

# Possible letters are those not in the incorrect set and earlier than 'S'
possible_letters = [letter for letter in string.ascii_uppercase if letter not in incorrect_letters and letter < 'S']

# The correct numbers are 1 and 8
correct_numbers = [1, 8]

# Print possible letters
print(possible_letters)