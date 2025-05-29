# Step 1: Division inside the parentheses
division_result = 7.512 / 7.152

# Step 2: Subtraction inside the parentheses
subtraction_inside_parentheses = 4.714 - division_result

# Step 3: Subtraction outside the parentheses
subtraction_outside_parentheses = subtraction_inside_parentheses - 9.191

# Step 4: Final division
final_result = subtraction_outside_parentheses / 3.225

# Print the final result rounded to 12 significant digits
print(f"{final_result:.12g}")