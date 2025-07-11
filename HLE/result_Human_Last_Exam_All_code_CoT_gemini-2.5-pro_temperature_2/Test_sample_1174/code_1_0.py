import sys

# The derived exponents from the physical analysis
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -2
n6 = -1.5

# The expression to be calculated is sum_{k=1 to 6} k*n_k
# This expands to: 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6

# Store coefficients and exponents in lists for clarity
k = [1, 2, 3, 4, 5, 6]
n = [n1, n2, n3, n4, n5, n6]

# Calculate the sum
total = 0
# Prepare a string to display the equation
equation_parts = []

for i in range(len(k)):
    term = k[i] * n[i]
    total += term
    # Use standard float formatting for consistency
    n_val_str = f"{n[i]:g}" if isinstance(n[i], float) and n[i] == int(n[i]) else f"{n[i]}"
    # Format the part of the equation string
    equation_parts.append(f"{k[i]}*({n_val_str})")

# On Python 3.8+ you might need to handle the output encoding for symbols explicitly
if sys.stdout.encoding != 'UTF-8':
    # A fallback for simple terminals
    equation_str = " + ".join(equation_parts)
else:
    # Use a proper summation symbol if supported
    equation_str = " + ".join(equation_parts)

# Print the full equation and the final result
# The format requested is to show the full equation leading to the result
print(f"Calculating the sum Σ k*n_k for k=1 to 6:")
print(f"{equation_str} = {total:g}")
