def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    # Store the digits to display them later
    digits = []
    for digit in str(n):
        s += int(digit)
        digits.append(digit)
    return s, digits

def find_next_in_sequence(start_after):
    """
    Finds the next number n > start_after where (n + sum of n's digits) is divisible by 6.
    """
    n = start_after + 1
    while True:
        s, digits = sum_digits(n)
        # The condition is that (n + S(n)) must be divisible by 6
        if (n + s) % 6 == 0:
            result = n + s
            # Print the equation as requested
            print(f"The next number is {n}.")
            print(f"The calculation is: {n} + ({' + '.join(digits)}) = {n} + {s} = {result}")
            print(f"This result, {result}, is divisible by 6.")
            return n
        n += 1

# The last number given in the sequence
last_number = 2352

# Find and print the next number
next_number = find_next_in_sequence(last_number)