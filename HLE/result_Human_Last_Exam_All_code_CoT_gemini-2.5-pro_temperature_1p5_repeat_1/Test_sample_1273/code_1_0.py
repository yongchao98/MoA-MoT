# Define the properties of the products
b2_price = 150
t1_price = 1

# Define the parameters of the optimal cutting strategy
num_b2_balls = 4
num_t1_cubes_from_waste_per_b2 = 12

# Calculate the total number of T1 cubes salvaged from the waste
total_t1_cubes = num_b2_balls * num_t1_cubes_from_waste_per_b2

# Calculate the total value from each type of product
value_from_b2 = num_b2_balls * b2_price
value_from_t1 = total_t1_cubes * t1_price

# Calculate the final maximum total value
total_value = value_from_b2 + value_from_t1

# Print the final equation showing the calculation
print(f"The optimal strategy is to cut {num_b2_balls} B2 balls and salvage {total_t1_cubes} T1 cubes from the waste material.")
print("The calculation for the maximum value is:")
print(f"({num_b2_balls} B2 balls * {b2_price}) + ({total_t1_cubes} T1 cubes * {t1_price}) = {value_from_b2} + {value_from_t1} = {total_value}")
