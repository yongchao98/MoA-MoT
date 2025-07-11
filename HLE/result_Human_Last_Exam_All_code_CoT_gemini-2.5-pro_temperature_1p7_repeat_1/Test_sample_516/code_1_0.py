import math

# The problem is to find the expected value of T = sum_{i=1 to m} Y_i / 3^i,
# where Y_i are i.i.d. U(0,1) and m is the smallest integer such that Y_m > Y_{m+1}.

# The derivation described above leads to the analytical result for the expected value:
# E[T] = 3 - 2 * e^(1/3)

# We can calculate this value using Python's math library.
# The numbers in the final equation are 3, 2, 1, and 3.

# Final equation for the expected value of T
val_3 = 3
val_2 = 2
power_numerator = 1
power_denominator = 3

expected_value = val_3 - val_2 * math.exp(power_numerator / power_denominator)

# Print the final equation and its numerical value.
# The format below shows each number in the final equation as requested.
print(f"The expected value E[T] is calculated by the equation: {val_3} - {val_2} * e^({power_numerator}/{power_denominator})")
print("E[T] =", expected_value)