# Set n, the size of the grid. You can change this value to see the result for a different grid size.
n = 4

# Explain the formula and the identified expressions for a and b.
print("To find the total number of squares on an n x n grid, we use the formula:")
print(f"Total Squares = Summation from m=1 to {n} of (n-m+1)^2 * m")
print("\nBy comparing this to the given expression Sum_{m=1 to n} a^2 * b, we can identify:")
print("a = n - m + 1")
print("b = m")

print(f"\nFor a grid of size n = {n}, we calculate the total number of squares by summing the terms for m from 1 to {n}:")

total_squares = 0
equation_terms = []
value_terms = []

# Loop through each possible bounding box size m.
for m in range(1, n + 1):
    # For each step of the summation, 'a' and 'b' are calculated based on n and the current m.
    a = n - m + 1
    b = m
    
    current_term_value = (a**2) * b
    total_squares += current_term_value
    
    # Build strings to display the final equation breakdown.
    # This shows each number used in the calculation.
    equation_terms.append(f"({a}^2 * {b})")
    value_terms.append(str(current_term_value))

# Present the final equation and its step-by-step evaluation.
final_equation = " + ".join(equation_terms)
final_values = " + ".join(value_terms)

print("\nThe final equation is:")
print(f"Total Squares = {final_equation}")
print(f"              = {final_values}")
print(f"              = {total_squares}")

print(f"\nThus, the total number of squares on a {n}x{n} grid is {total_squares}.")
