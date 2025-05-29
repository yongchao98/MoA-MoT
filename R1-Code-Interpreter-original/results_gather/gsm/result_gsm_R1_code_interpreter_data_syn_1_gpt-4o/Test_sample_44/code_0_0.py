# Initial number of shoes
initial_shoes = 3495007

# Shoes gained on Monday
monday_gain = 5 * 2

# Shoes gained on Wednesday
wednesday_gain = 15 * 2

# Shoes gained on Friday
friday_gain = 30 * 2

# Shoes lost on Saturday
saturday_loss = 180

# Total shoes on Sunday
total_shoes = initial_shoes + monday_gain + wednesday_gain + friday_gain - saturday_loss

print(total_shoes)