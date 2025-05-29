# Possible numbers and letters
numbers = ['2', '7']
letters = ['J', 'M']

# Possible positions for numbers and letters
# From guess 11, 7 and 2 are correct but in the wrong positions, so the correct order is 27.
number_positions = ['2', '7']

# From guess 3, one letter is correct and in the correct position.
# From guess 5, one letter is correct and in the correct position.
# The correct letters are J and M, and they must be in the same position as S or X in another guess.
letter_positions = ['J', 'M']

# Final password
password = number_positions + letter_positions
print(password)