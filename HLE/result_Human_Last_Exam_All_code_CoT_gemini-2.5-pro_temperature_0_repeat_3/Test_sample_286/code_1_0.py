import math

# The core of the problem reduces to the inequality:
# 1 - 10^(-M) >= 1/2
# which simplifies to M >= log10(2)

# Let's define the numbers in the final equation
# 2 * B >= 1  => 2 * (1 - 10**(-M)) >= 1
val_2 = 2
val_1 = 1
val_10 = 10

print(f"The length of the interval for x_i is 2*B, where B = 1 - 10^(-M).")
print(f"A single step 'a_i' can be at most the length of the interval, so a_i <= 2*B.")
print(f"Since a_i can be 1, we must have 1 <= 2*B.")
print(f"This gives the inequality: {val_1} <= {val_2} * ({val_1} - {val_10}**(-M))")

# This simplifies to 1/2 >= 10**(-M), or M >= log10(2)
log_base = 10
argument = 2
lower_bound = math.log10(argument)

print(f"\nSolving for M, we get M >= log_{log_base}({argument})")
print(f"The value of log_{log_base}({argument}) is approximately {lower_bound}.")

# M must be the smallest positive integer satisfying this condition.
M = math.ceil(lower_bound)

print(f"The smallest positive integer M satisfying M >= {lower_bound} is {M}.")
