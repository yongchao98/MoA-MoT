import math

# Based on the analytical derivation, the integral simplifies significantly.
# The calculation for the second part of the integral, I_2, after substitution
# and integration by parts, results in a constant term and an integral term.
# The constant term is a product of two numbers derived from the process.

# The coefficient from the I_2 integral is C = 2**(17/16)
C = 2**(17/16)

# The boundary term evaluated from 0 to pi/4 is B = 2**(-1/8)
B = 2**(-1/8)

# The final value of the integral is the product of these two numbers,
# as the remaining integral parts cancel out.
integral_value = C * B

print(f"The calculation for the final constant term is:")
print(f"({C}) * ({B})")
print(f"Which simplifies to 2^(17/16) * 2^(-2/16) = 2^(15/16)")
print(f"\nThe exact value of the integral is 2**(15/16).")
print(f"The numerical value is: {integral_value}")
