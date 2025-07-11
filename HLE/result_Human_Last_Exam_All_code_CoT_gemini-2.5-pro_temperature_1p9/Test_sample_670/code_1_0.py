import sys
# This is a symbolic representation, so we use strings.
# The `μ` (mu) operator represents the least fixed-point for a recursive definition, used for loops.
# It binds a variable, which we call `a` as requested.
# The `φ` (phi) function represents the merging of values at the start of the loop.
# The first argument to `φ` is the initial value from before the loop.
# The second argument is the value from the loop's back-edge (the previous iteration).
initial_value = 1
# Let 'a' be the bound variable representing the value of the variable 'a' from the previous iteration.
# Let 'i' be the loop's induction variable.
# The expression for the value computed in the loop body is "a * i".
print(f"μ a . φ({initial_value}, a * i)")