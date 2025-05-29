# Initialize the list of fruits with None
fruits = [None] * 7

# Place the known positions
fruits[1] = 'watermelons'  # Second-most expensive
fruits[3] = 'oranges'      # Fourth-most expensive
fruits[5] = 'pears'        # Second-cheapest

# Determine the positions of the remaining fruits
# We know:
# - Peaches are more expensive than cantaloupes and less expensive than mangoes
# - Cantaloupes are more expensive than kiwis

# Let's try to place the remaining fruits
# Possible positions for mangoes, peaches, cantaloupes, kiwis: 0, 2, 4, 6

# Since watermelons are in position 1, mangoes must be in position 0 (most expensive)
fruits[0] = 'mangoes'

# Peaches must be less expensive than mangoes and more expensive than cantaloupes
# So, peaches can be in position 2
fruits[2] = 'peaches'

# Cantaloupes are more expensive than kiwis, so cantaloupes can be in position 4
fruits[4] = 'cantaloupes'

# Finally, kiwis must be the cheapest
fruits[6] = 'kiwis'

# Print the order of fruits
print(fruits)