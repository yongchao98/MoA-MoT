# Plan: Solve the word puzzle by extracting numbers from the text and summing them up.

# 1. Extract the numbers from the text.
# The year is 2008. We will use its individual digits.
year = 2008
year_digits = [int(digit) for digit in str(year)] # This gives [2, 0, 0, 8]

# "single-turn" implies the number 1.
turn_count = 1

# "fifth position" implies the number 5.
position_number = 5

# "Act l" implies the number 1.
act_number = 1

# 2. Sum all the extracted numbers.
total_sum = sum(year_digits) + turn_count + position_number + act_number

# 3. Print the full equation showing each number.
# We build a list of all the numbers we are adding.
all_numbers = year_digits + [turn_count, position_number, act_number]

# We create the string for the equation.
equation_str = " + ".join(map(str, all_numbers))

# 4. Print the final result.
print(f"The final equation is: {equation_str} = {total_sum}")