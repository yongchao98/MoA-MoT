# The pattern from the Bansenshukai scroll: ⬤○○⬤⬤⬤⬤○⬤⬤⬤⬤⬤
pattern_string = "BOOBBBBOBBBBB"

# Count the number of black circles (represented by 'B')
num_black_circles = pattern_string.count('B')

# Count the number of white circles (represented by 'O')
num_white_circles = pattern_string.count('O')

# Calculate the total number of circles
total_circles = num_black_circles + num_white_circles

# Print the final equation showing the counts
print(f"Counting the symbols reveals an equation:")
print(f"{num_black_circles} + {num_white_circles} = {total_circles}")