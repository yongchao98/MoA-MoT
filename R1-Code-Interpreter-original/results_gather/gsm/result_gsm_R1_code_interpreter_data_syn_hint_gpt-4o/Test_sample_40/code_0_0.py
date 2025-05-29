# Total number of hats
total_hats = 716646

# Calculate the number of red hats
red_hats = total_hats / 2

# Calculate the remaining hats after removing red hats
remaining_hats = total_hats - red_hats

# Calculate the number of blue hats
blue_hats = (3/5) * remaining_hats

# Calculate the number of white hats
white_hats = remaining_hats - blue_hats

# Output the number of white hats
print(int(white_hats))