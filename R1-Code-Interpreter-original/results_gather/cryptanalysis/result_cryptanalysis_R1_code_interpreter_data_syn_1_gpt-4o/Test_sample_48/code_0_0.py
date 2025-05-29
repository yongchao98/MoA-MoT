# Possible numbers and letters based on deductions
possible_numbers = [6, 9]
possible_letters = ['X', 'Z']

# Since we know the numbers and letters, we need to find the correct positions
# Based on the feedback, we know:
# - 6 and 9 are in the wrong positions in the guesses
# - Z is in the correct position in guess 2
# - X is in the wrong position in guess 5

# The correct combination should be:
# - 9 in the first position (since it was in the wrong position in guess 6)
# - 6 in the second position (since it was in the wrong position in guesses 2 and 4)
# - Z in the third position (since it was correct in guess 2)
# - X in the fourth position (since it was in the wrong position in guess 5)

combination = [9, 6, 'Z', 'X']
print(combination)