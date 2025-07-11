import decimal

# 1. Set the precision for decimal calculations.
# The number 10^100 has 101 digits. To get an accurate remainder
# when dividing by pi, the precision must be greater than 101.
# We set it to 125 for a safe margin.
decimal.getcontext().prec = 125

# 2. Define a high-precision value of pi as a string.
# This string contains pi to over 120 decimal places.
pi_str = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446"
pi = decimal.Decimal(pi_str)

# 3. Define the number x = 10^100 as a Decimal object.
x = decimal.Decimal(10)**100

# 4. Use the periodicity of tan: tan(x) = tan(x mod pi).
# We calculate the argument r = x mod pi.
argument = x % pi

# 5. Calculate the tangent of the small remainder 'argument'.
result = argument.tan()

# 6. Output the numbers involved in the final equation.
print(f"The calculation is based on the identity: tan(x) = tan(x mod pi)")
print(f"The argument (10^100 mod pi) is calculated to be: {argument}")
print(f"The final result for tan(10^100) is therefore tan({argument}), which equals: {result}")

# 7. Extract and print the first 3 digits after the comma.
# To get the digits, we can take the absolute value of the result,
# multiply by 1000, convert to an integer, and take the result modulo 1000.
first_three_digits = int(abs(result) * 1000) % 1000

print(f"\nThe first 3 digits after the comma are: {first_three_digits}")