# Values of the exponents derived from the physical analysis
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -1
n6 = -0.5

# The coefficients for the sum
k = [1, 2, 3, 4, 5, 6]
n = [n1, n2, n3, n4, n5, n6]

# Calculate the sum of k*n_k
total_sum = 0
equation_parts = []
for i in range(len(k)):
    term = k[i] * n[i]
    total_sum += term
    # Format each term for the equation string
    # Show n_k as a number, handling floats gracefully
    n_val_str = f"{n[i]:g}" 
    equation_parts.append(f"{k[i]}*({n_val_str})")

# Join the parts to form the full equation string
equation_str = " + ".join(equation_parts)

# Print the full equation and the final result
print(f"The calculation is: Σ k*n_k = {equation_str}")
print(f"The final value is: {total_sum:g}")
