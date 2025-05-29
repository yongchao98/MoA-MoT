# Original price of the radiator
original_price = 1199484

# Discount rate
discount_rate = 0.80

# Calculate the discounted price
discounted_price = original_price * (1 - discount_rate)

# Mechanic's hourly rate and hours worked
hourly_rate = 50
hours_worked = 3

# Calculate the installation cost
installation_cost = hourly_rate * hours_worked

# Calculate the total cost
total_cost = discounted_price + installation_cost

# Round the total cost to two decimal places
total_cost_rounded = round(total_cost, 2)

# Print the rounded total cost
print(total_cost_rounded)