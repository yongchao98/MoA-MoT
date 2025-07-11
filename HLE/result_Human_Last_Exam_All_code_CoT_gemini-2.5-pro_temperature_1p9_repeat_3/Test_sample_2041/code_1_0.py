# The task is to find the number of "shallow" polymorphic functions, which boils down to
# counting the number of boolean functions of 4 variables.

# Let k be the number of variables. Based on the analysis, we have 4 independent boolean
# values derived from applying p, so k=4.
k = 4

# The final result is calculated by the equation 2^(2^k).
# We are asked to output each number in this final equation.
outer_base = 2
inner_base = 2

# First, calculate the size of the truth table for a function of k variables.
truth_table_size = inner_base**k

# Then, calculate the total number of functions. This is 2 to the power of the
# truth table size, since each row can result in True or False.
total_functions = outer_base**truth_table_size

print(f"The number of variables (k) is {k}.")
print(f"The number of functions is given by the equation: {outer_base} ** ({inner_base} ** {k})")
print(f"The calculation proceeds in two steps:")
print(f"1. Inner exponentiation: {inner_base} ** {k} = {truth_table_size}")
print(f"2. Outer exponentiation: {outer_base} ** {truth_table_size} = {total_functions}")
print(f"Thus, there are {total_functions} extensionally distinct functions.")