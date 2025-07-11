import sys

# Determined exponents from the physics of thermal magnetic noise
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -2
n6 = -0.5

# The coefficients are the 'k' values from the summation
coeffs = [1, 2, 3, 4, 5, 6]
exponents = [n1, n2, n3, n4, n5, n6]

# Calculate the sum
total = sum(k * n for k, n in zip(coeffs, exponents))

# Output the equation and the result
# The problem asks to output the equation with each number.
# We build the string for the equation first.
equation_parts = []
for k, n in zip(coeffs, exponents):
    # Handle negative numbers for nice printing
    if n < 0:
        equation_parts.append(f"{k}*({n})")
    else:
        equation_parts.append(f"{k}*{n}")

equation_str = " + ".join(equation_parts)
print(f"The calculation is: {equation_str}")
print(f"The result is: {total}")

# The required answer format.
# Note that we are using `sys.stdout.write` to avoid the default newline from `print`
# to match the desired format `<<<answer content>>>` precisely.
sys.stdout.write(f'<<<{total}>>>')