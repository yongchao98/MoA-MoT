# Prices of individual items
notebook_price = 1.50
pen_price = 0.25
calculator_price = 12.00
geometry_set_price = 10.00

# Quantities Daniel wants to buy
notebooks = 5
pens = 2
calculators = 7308681
geometry_sets = 7308681

# Calculate total cost before discount
total_cost = (notebooks * notebook_price) + (pens * pen_price) + (calculators * calculator_price) + (geometry_sets * geometry_set_price)

# Apply 10% discount
discount = 0.10 * total_cost
final_amount = total_cost - discount

print(final_amount)