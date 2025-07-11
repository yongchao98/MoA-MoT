# This script calculates the number of extensionally distinct functions as derived
# from the step-by-step analysis.
# The problem reduces to finding the number of boolean functions of 4 variables.

# The number of independent atomic propositions that can be formed.
num_atoms = 4

# The number of possible truth-value assignments for these atoms.
num_truth_assignments = 2**num_atoms

# The total number of distinct functions is the number of possible boolean
# functions of `num_atoms` variables, which is 2 to the power of the
# number of possible truth assignments.
num_functions = 2**num_truth_assignments

# The final equation is 2**(2**4) = 65536.
# As requested, the following prints each number in that final equation.
print("The final calculation is based on the equation: 2**(2**4) = 65536")
print("The numbers that form this equation are:")
print(2)
print(2)
print(4)
print(num_functions)
