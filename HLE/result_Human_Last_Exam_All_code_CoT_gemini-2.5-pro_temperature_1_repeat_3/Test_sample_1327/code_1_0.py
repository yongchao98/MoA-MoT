# The last term in the given sequence.
last_term = 119

# The last number that was subtracted in the pattern.
# For 119, the calculation was 3 * 51 - 34 = 119. So, 34 is the last subtrahend.
last_subtrahend = 34

# The difference between the last two subtrahends (18 and 34) was 16.
last_difference = 16

# The pattern of the differences is that they double each time.
next_difference = last_difference * 2

# Calculate the next number to subtract.
next_subtrahend = last_subtrahend + next_difference

# Calculate the next term in the main sequence.
next_term = 3 * last_term - next_subtrahend

# Print the final equation with each number.
print(f"The calculation for the next term is:")
print(f"3 * {last_term} - {next_subtrahend} = {next_term}")