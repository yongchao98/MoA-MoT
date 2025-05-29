# Prices of individual items
price_notebook = 1.50
price_pen = 0.25
price_calculator = 12.00
price_geometry_set = 10.00

# Quantities Daniel wants to buy
quantity_notebook = 5
quantity_pen = 2
quantity_calculator = 7308681
quantity_geometry_set = 7308681

# Calculate total cost for each item
total_cost_notebooks = price_notebook * quantity_notebook
total_cost_pens = price_pen * quantity_pen
total_cost_calculators = price_calculator * quantity_calculator
total_cost_geometry_sets = price_geometry_set * quantity_geometry_set

# Calculate total cost before discount
total_cost_before_discount = (total_cost_notebooks + total_cost_pens +
                              total_cost_calculators + total_cost_geometry_sets)

# Apply 10% discount
discount = 0.10 * total_cost_before_discount
total_cost_after_discount = total_cost_before_discount - discount

# Print the final amount Daniel has to spend
print(total_cost_after_discount)