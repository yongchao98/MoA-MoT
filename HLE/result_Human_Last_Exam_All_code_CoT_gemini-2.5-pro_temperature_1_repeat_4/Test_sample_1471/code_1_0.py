def convert_base17_to_decimal(base17_digits):
    """
    Converts a number represented as a list of base-17 digits into its decimal (base-10) equivalent.
    """
    decimal_value = 0
    power = len(base17_digits) - 1
    for digit in base17_digits:
        if not 0 <= digit < 17:
            raise ValueError("Invalid digit for base-17 representation.")
        decimal_value += digit * (17 ** power)
        power -= 1
    return decimal_value

# The problem uses 'септендецимал', which means base-17.
# The symbols 𒌋 and 𒁹 correspond to the values 10 and 1.

# First number (𒌋𒌋𒌋) is represented by digits [10, 10, 10] in base-17.
num1_digits = [10, 10, 10]
num1_decimal = convert_base17_to_decimal(num1_digits)

# Second number (𒌋𒌋𒁹𒁹𒁹) is represented by digits [10, 10, 1, 1, 1] in base-17.
num2_digits = [10, 10, 1, 1, 1]
num2_decimal = convert_base17_to_decimal(num2_digits)

# Perform the subtraction
result = num1_decimal - num2_decimal

# Print the final equation using modern numbers
print(f"The equation in modern numbers is:")
print(f"{num1_decimal} - {num2_decimal} = {result}")
