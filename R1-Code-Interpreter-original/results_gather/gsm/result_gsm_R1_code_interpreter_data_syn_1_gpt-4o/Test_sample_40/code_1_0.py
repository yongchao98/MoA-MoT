# Total number of hats
total_hats = 716646

# Calculate the number of red hats using integer division
red_hats = total_hats // 2

# Calculate the remaining hats after removing red hats
remaining_hats = total_hats - red_hats

# Calculate the number of blue hats using integer division
blue_hats = (3 * remaining_hats) // 5

# Calculate the number of white hats
white_hats = remaining_hats - blue_hats

# Output the number of white hats
print(white_hats)