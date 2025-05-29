# Prices of individual items
price_notebook = 1.50
price_pen = 0.25
price_calculator = 12.00
price_geometry_set = 10.00

# Quantities Daniel wants to buy
quantity_notebooks = 5
quantity_pens = 2
quantity_calculators = 7308681
quantity_geometry_sets = 7308681

# Calculate total cost before discount
total_cost = (quantity_notebooks * price_notebook) + \
             (quantity_pens * price_pen) + \
             (quantity_calculators * price_calculator) + \
             (quantity_geometry_sets * price_geometry_set)

# Apply 10% discount
discounted_total_cost = total_cost * (1 - 0.10)

# Output the final amount
print(discounted_total_cost)