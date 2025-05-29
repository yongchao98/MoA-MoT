# Check if the position (1, 4) fits into the pattern
def check_position(row, col):
    difference = row - col
    return difference in differences

# Check the position of the '❔'
position_check = check_position(1, 4)

# Print the result of the check
print(position_check)