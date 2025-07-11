import numpy as np

# Let H(x,t) be the result of the fractional differentiation.
# From the analytical derivation, we found:
# H(x,t) = -(exp(-(x + 6*t)) + exp(-2*(x + 6*t)))

# The problem asks to evaluate this quantity at x = 6t.
# This results in H(6t, t) = -(exp(-12*t) + exp(-24*t)).
# To get a single numerical value, we evaluate at the origin (x=0, t=0),
# which lies on the line x = 6t.

x = 0
t = 0

# The expression to be calculated is -(e^-(x+6t) + e^-2(x+6t))
# We substitute the values of x and t.

val1 = -(x + 6 * t)
val2 = -2 * (x + 6 * t)

# Calculate the final result
result = -(np.exp(val1) + np.exp(val2))

# Print the final equation with all numbers
print(f"Result = -(exp(-({x} + 6*{t})) + exp(-2*({x} + 6*{t})))")
print(f"Result = -(exp({val1}) + exp({val2}))")
print(f"Result = -({np.exp(val1)} + {np.exp(val2)})")
print(f"Result = {result}")

# Final Answer
# The final answer is wrapped in <<<>>>
print("<<<{}>>>".format(result))