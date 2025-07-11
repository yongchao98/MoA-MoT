# Set the value of n for the n x n grid
n = 5

print(f"Calculating the number of squares for an n x n grid where n = {n}")
print("The formula is sum_{m=1 to n} a^2 * b")
print("Based on our derivation, this becomes sum_{m=1 to n} (n-m+1)^2 * m")
print(f"where a = n-m+1 and b = m.")

total_squares = 0
equation_parts = []
numeric_parts = []

# Loop through the summation from m = 1 to n
for m in range(1, n + 1):
    # For each m, a = n-m+1 and b = m
    a = n - m + 1
    b = m
    
    # The term in the sum is a^2 * b
    term_value = (a**2) * b
    total_squares += term_value
    
    # Store the string representation of each part of the equation
    equation_parts.append(f"{a}^2*{b}")
    numeric_parts.append(str(term_value))

# Join the parts to form the full equation string
equation_str = " + ".join(equation_parts)
numeric_str = " + ".join(numeric_parts)

print("\nThe full equation with the numbers filled in is:")
print(f"Total Squares = {equation_str}")
print(f"              = {numeric_str}")
print(f"              = {total_squares}")

print("\nTherefore, to make the expression equal to the number of squares, 'a' and 'b' must be:")
print("a = n - m + 1")
print("b = m")