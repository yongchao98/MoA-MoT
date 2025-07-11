# Step 1: Assign the values of the exponents based on the physical analysis.
n1 = 3
n2 = 1
n3 = 2
n4 = 2
n5 = -2
n6 = -3/2

# Step 2: Define the coefficients for the sum.
k = [1, 2, 3, 4, 5, 6]
n = [n1, n2, n3, n4, n5, n6]

# Step 3: Calculate the sum and print the equation.
total_sum = 0
equation_parts = []
for i in range(6):
    term = k[i] * n[i]
    total_sum += term
    # Format each term for the equation string
    # Use int() to display integer results without .0
    n_val_str = str(int(n[i])) if n[i] == int(n[i]) else str(n[i])
    equation_parts.append(f"{k[i]}*({n_val_str})")

# Step 4: Print the full equation and the final result.
equation_str = " + ".join(equation_parts)
print(f"The calculation is: {equation_str}")
print(f"The final sum is: {int(total_sum)}")
