# Prices of each item
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
total_cost_before_discount = (notebooks * notebook_price) + (pens * pen_price) + (calculators * calculator_price) + (geometry_sets * geometry_set_price)

# Apply 10% discount
total_cost_after_discount = total_cost_before_discount * 0.90

# Output the final amount
print(total_cost_after_discount)