def sum_digits(n):
    """Calculates the sum of the digits of a number."""
    s = 0
    # Store the original number to print later
    original_n = n
    digits = []
    while n > 0:
        digit = n % 10
        s += digit
        digits.append(str(digit))
        n //= 10
    
    # For printing purposes, let's create the equation string
    # E.g., for 123, it returns (6, "1 + 2 + 3 = 6")
    equation = " + ".join(reversed(digits)) + f" = {s}"
    return s, equation

def is_harshad(n):
    """Checks if a number is a Harshad number."""
    if n == 0:
        return False, None
    s, _ = sum_digits(n)
    return n % s == 0, s

def find_next_number_in_sequence(start_num):
    """
    Finds the next number in the sequence starting from start_num + 1.
    The sequence is Harshad numbers whose sum of digits is also a Harshad number.
    """
    num = start_num + 1
    while True:
        # Check if num is a Harshad number
        is_num_harshad, num_digit_sum = is_harshad(num)
        
        if is_num_harshad:
            # If it is, check if its sum of digits is also a Harshad number
            is_sum_harshad, _ = is_harshad(num_digit_sum)
            if is_sum_harshad:
                # We found the number
                return num
        num += 1

# The last number in the provided sequence
last_number = 2352

# Find the next number
next_number = find_next_number_in_sequence(last_number)

print(f"The last number in the sequence is {last_number}.")
print("The rule is that the number must be a Harshad number, and its sum of digits must also be a Harshad number.")
print(f"Searching for the next number after {last_number}...")
print(f"\nFound the next number in the sequence: {next_number}\n")
print(f"Verification for {next_number}:")

# Verification Step 1: Check if next_number is a Harshad number
s_num, eq_num = sum_digits(next_number)
print(f"1. Check if {next_number} is a Harshad number (divisible by the sum of its digits).")
print(f"   - Sum of digits for {next_number}: {eq_num}")
print(f"   - Division check: {next_number} / {s_num} = {next_number // s_num}")
print(f"   - Since the division has no remainder, {next_number} is a Harshad number.")

# Verification Step 2: Check if the sum of digits is a Harshad number
s_sum, eq_sum = sum_digits(s_num)
print(f"\n2. Check if the sum of its digits ({s_num}) is also a Harshad number.")
print(f"   - Sum of digits for {s_num}: {eq_sum}")
print(f"   - Division check: {s_num} / {s_sum} = {s_num // s_sum}")
print(f"   - Since the division has no remainder, {s_num} is also a Harshad number.")

print(f"\nBoth conditions are met. Therefore, the next number in the sequence is {next_number}.")

print(f"\n<<<2358>>>")