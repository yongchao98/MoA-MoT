import math

# The problem simplifies to finding the integer part of 10^4 * L,
# where L is the limit of ln(s_n)/n.
# We found that L = ln(R), where R is the dominant root of the
# characteristic equation of the recurrence for a_n.
# The dominant root is R = 2 + sqrt(7).

# The value to calculate is 10^4 * ln(2 + sqrt(7)).
R = 2 + math.sqrt(7)
limit_val = math.log(R)
final_result = 10000 * limit_val

# We need the integer part of this result.
integer_part = math.floor(final_result)

# As requested, printing the number from the final equation
# The "equation" here is floor(10^4 * ln(2 + sqrt(7)))
# where each number needs to be printed.
print(f"10000")
print(f"2")
print(f"7")
print(f"The final integer value is: {integer_part}")