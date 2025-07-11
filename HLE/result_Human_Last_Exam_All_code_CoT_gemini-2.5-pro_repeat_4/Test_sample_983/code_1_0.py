def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    while n:
        s += n % 10
        n //= 10
    return s

# The last number in the sequence
last_number = 2352

# Calculate the sum of its digits
s_digits = sum_digits(last_number)

# The rule to find the increment is 2 * sum_of_digits
increment = 2 * s_digits

# Calculate the next number
next_number = last_number + increment

# Print the equation
print(f"{last_number} + 2 * {s_digits} = {next_number}")
print(f"{last_number} + {increment} = {next_number}")

# The final answer
# print(f"<<<{next_number}>>>")