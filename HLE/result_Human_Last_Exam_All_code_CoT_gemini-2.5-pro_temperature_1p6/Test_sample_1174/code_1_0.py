# Values of the exponents n_k based on physics derivation
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -2
n6 = -0.5

# The list of coefficients k from 1 to 6
k_values = [1, 2, 3, 4, 5, 6]

# The list of exponents n_k
n_values = [n1, n2, n3, n4, n5, n6]

# Calculate the sum of k * n_k
total_sum = 0
equation_parts = []
for k, n in zip(k_values, n_values):
    term = k * n
    total_sum += term
    # Format the number to be an integer if it's a whole number
    if n == int(n):
        n_str = str(int(n))
    else:
        n_str = str(n)
        
    equation_parts.append(f"{k}*({n_str})")

# Print the full equation
equation_str = " + ".join(equation_parts)
print(f"The calculation is: {equation_str}")

# Print the final result
# Format the number to be an integer if it's a whole number
if total_sum == int(total_sum):
    result_str = str(int(total_sum))
else:
    result_str = str(total_sum)
print(f"The value of the sum is: {result_str}")
