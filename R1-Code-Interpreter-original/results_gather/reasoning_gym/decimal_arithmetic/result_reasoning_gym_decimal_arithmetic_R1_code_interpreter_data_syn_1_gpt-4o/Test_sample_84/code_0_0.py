# Step 1: Calculate the innermost division
inner_division = 4.491 / 9.387

# Step 2: Divide 7.877 by the result of the inner division
division_result = 7.877 / inner_division

# Step 3: Calculate the multiplication
multiplication_result = 1.047 * 8.131

# Step 4: Add the results from steps 2 and 3
inner_addition = division_result + multiplication_result

# Step 5: Add 0.144 to the result from step 4
final_result = 0.144 + inner_addition

# Print the final result rounded to 12 significant digits
print(f"{final_result:.12g}")