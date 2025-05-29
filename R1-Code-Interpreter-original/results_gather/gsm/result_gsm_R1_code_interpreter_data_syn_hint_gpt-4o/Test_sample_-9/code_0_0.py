# Initial number of shoes
initial_shoes = 3495007

# Shoes gained
monday_gain = 5 * 2
wednesday_gain = 15 * 2
friday_gain = 30 * 2

# Total shoes gained
total_gain = monday_gain + wednesday_gain + friday_gain

# Shoes lost
saturday_loss = 180

# Net change in shoes
net_change = total_gain - saturday_loss

# Final number of shoes
final_shoes = initial_shoes + net_change

print(final_shoes)