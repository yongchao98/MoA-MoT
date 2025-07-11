# This script calculates the number of extensionally distinct functions
# induced by "shallow" expressions in the given lambda calculus setting.
# The final result is derived by calculating 2^(2^(k^m)), where:
# k is the number of available base predicates.
# m is the number of arguments p takes.

# Step 1: Determine the number of valid base predicates (k).
# An argument to `p` must be a `p-free` predicate of type X->Bool.
# In a standard simply typed lambda calculus with an opaque type X,
# there are only two such extensionally distinct predicates:
# 1. The predicate that always returns True (`lambda x: True`).
# 2. The predicate that always returns False (`lambda x: False`).
k = 2

# Step 2: Determine the arity of `p` (m).
# The type of `p` is PPPX, meaning p(predicate)(predicate)(predicate), so it takes 3 arguments.
m = 3

# Step 3: Determine the number of fundamental boolean signals (n = k^m).
# A shallow expression can "probe" `p` by applying it to any combination of the k base predicates.
# The number of such combinations (probes) is k^m.
# Each probe, like `p(P_True, P_False, P_True)(x)`, yields a single boolean value.
# These values are the fundamental signals.
n = k**m

# Step 4: Calculate the total number of functions.
# A shallow expression `e` can represent any boolean function of these `n` signals.
# The number of boolean functions of n variables is 2^(2^n).
# The final result is 2^(2^n) = 2^(2^(2^3)).

# We calculate the final result and print the numbers involved in the equation.
power_of_2_n = 2**n
final_result = 2**power_of_2_n

print(f"The number of base p-free predicates (k) is: {k}")
print(f"The number of arguments `p` takes (m) is: {m}")
print(f"The number of distinct boolean signals (n = k^m) is {k}^{m}, which is: {n}")
print(f"The total number of functions is 2^(2^n). We calculate 2^{n}, which is: {power_of_2_n}")
print(f"Finally, the total number of distinct functions is 2^{power_of_2_n}.")
print("\nFinal Answer:")
print(final_result)