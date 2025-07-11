import math

def catalan(n):
    """Calculates the n-th Catalan number."""
    if n < 0:
        return 0
    # C_n = comb(2n, n) // (n + 1)
    return math.comb(2 * n, n) // (n + 1)

def calculate_dn(n):
    """Calculates the value of d_n based on the given formula."""
    if n == 1:
        return 2
    elif n > 1:
        return catalan(n)
    else:
        return 0

# Calculate d_n for n from 1 to 8
d_values = [calculate_dn(n) for n in range(1, 9)]

# Calculate the product
product = 1
for d in d_values:
    product *= d

# Format the output string
equation_str = " * ".join(map(str, d_values))

# Print the final equation and the result
print(f"The values for d_n from n=1 to 8 are: {d_values}")
print("The product is:")
print(f"{equation_str} = {product}")