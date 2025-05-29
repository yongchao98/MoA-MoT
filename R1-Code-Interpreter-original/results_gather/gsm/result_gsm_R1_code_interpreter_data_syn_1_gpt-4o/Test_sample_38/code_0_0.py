# Given data
sauce_volume = 32  # ounces
can_volume = 16  # ounces per can
tomatoes_per_can = 3  # tomatoes per can

# Calculate the original volume of tomatoes before cooking
original_volume = sauce_volume * 2

# Calculate the number of cans needed
number_of_cans = original_volume / can_volume

# Calculate the number of tomatoes used
number_of_tomatoes = number_of_cans * tomatoes_per_can

print(number_of_tomatoes)