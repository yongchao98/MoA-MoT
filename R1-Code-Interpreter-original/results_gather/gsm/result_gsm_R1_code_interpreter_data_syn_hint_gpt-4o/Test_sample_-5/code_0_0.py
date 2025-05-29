# Original cost of the radiator
original_cost = 1199484

# Discount percentage
discount_percentage = 80

# Cost after discount
cost_after_discount = original_cost * (1 - discount_percentage / 100)

# Mechanic's hourly rate
hourly_rate = 50

# Hours taken for installation
hours = 3

# Installation cost
installation_cost = hourly_rate * hours

# Total cost
total_cost = cost_after_discount + installation_cost

print(total_cost)