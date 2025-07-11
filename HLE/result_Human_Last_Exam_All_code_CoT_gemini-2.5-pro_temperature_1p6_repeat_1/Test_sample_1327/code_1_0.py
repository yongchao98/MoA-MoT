# The given sequence
sequence = [2, 11, 23, 51, 119]

# The last known number in the sequence
last_number = sequence[-1]

# The sequence of subtracted values is 10, 18, 34.
# The differences between these values are 8 (18-10) and 16 (34-18).
# The next difference is 16 * 2 = 32.
# So, the next value to subtract is 34 + 32 = 66.
k = 66

# The pattern is: next_number = 3 * previous_number - k
multiplier = 3
next_number = multiplier * last_number - k

# Output the equation and the result
print(f"The pattern is next_number = 3 * previous_number - k.")
print(f"The next k is calculated as 34 + (16 * 2) = {k}")
print(f"So, the equation is: {multiplier} * {last_number} - {k} = {next_number}")