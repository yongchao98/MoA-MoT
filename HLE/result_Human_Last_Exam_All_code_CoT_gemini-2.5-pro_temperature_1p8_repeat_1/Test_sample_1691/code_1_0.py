import math

# Parameters from the integral's denominator.
# The term with the lowest power of x is 9.0 * x^5.0.
# p is the coefficient and n is the power of this dominant term.
p = 9.0
n = 5.0

# The approximation for the integral has the form: I(epsilon) approx C * epsilon**(-a)

# The exponent 'a' is determined by the dominant power 'n'.
a = (n - 1.0) / n

# The coefficient 'C' is determined by 'p' and 'n'.
# The formula for C involves a standard definite integral, whose value is (pi/n) / sin(pi/n).
integral_value = (math.pi / n) / math.sin(math.pi / n)
C = p**(-1.0 / n) * integral_value

# The script will now print the derivation and the final formula clearly.

print("The analytical formula for I(epsilon) for small epsilon is of the form:")
print("I(epsilon) approx C * epsilon**(-a)")
print("\nThe values for the constants 'a' and 'C' are calculated as follows:")

# Explain and print the value of 'a'.
print(f"\nThe exponent 'a' is derived from the lowest power of x (n={n}):")
print(f"a = (n - 1) / n = ({n} - 1) / {n} = {a}")

# Explain and print the value of 'C'.
print(f"\nThe coefficient 'C' is derived from the lowest power of x (n={n}) and its coefficient (p={p}):")
print(f"C = p**(-1/n) * ( (pi/n) / sin(pi/n) )")
print(f"C = {p}**(-1/{n}) * ( (pi/{n}) / sin(pi/{n}) )")
# Print the intermediate numerical values for clarity
print(f"C approx {p**(-1.0/n):.4f} * {integral_value:.4f} = {C:.4f}")

# Print the final numerical formula as requested, with all numbers included.
print("\nThus, the final approximate analytical formula is:")
print(f"I(epsilon) approx {C:.4f} * epsilon**(-{a})")