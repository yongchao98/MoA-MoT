# Given number of green eggs
green_eggs = 7248186

# Calculate the number of blue eggs
blue_eggs = 2 * green_eggs

# Calculate the number of pink eggs
pink_eggs = blue_eggs - 1

# Calculate the number of yellow eggs
yellow_eggs = pink_eggs // 3  # Using integer division

# Calculate the total number of eggs
total_eggs = green_eggs + blue_eggs + pink_eggs + yellow_eggs

# Output the total number of eggs
print(total_eggs)