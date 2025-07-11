import math

# The value of t is ln(10^20)
# For this t, e^t = 10^20 and e^-t = 10^-20

# We can define these values directly
et = 1e20
e_neg_t = 1e-20

# The analytical solution to the boundary-value problem is:
# x(t) = (2 - e^(-t)) / cosh(t)
# which can be rewritten in terms of exponentials as:
# x(t) = 2 * (2 - e^(-t)) / (e^t + e^(-t))

# Calculate the numerator and denominator using these values.
# Note: Standard double-precision floating-point arithmetic has limited precision.
# The term 'e_neg_t' is so small compared to '2' and 'et' that it will be rounded off
# in the addition/subtraction operations.
# 2 - 1e-20 will be evaluated as 2.0
# 1e20 + 1e-20 will be evaluated as 1e20
numerator = 2 * (2 - e_neg_t)
denominator = et + e_neg_t

# Calculate the final result
result = numerator / denominator

print("The solution to the boundary-value problem is x(t) = 2 * (2 - exp(-t)) / (exp(t) + exp(-t)).")
print("We need to evaluate x(t) for t = ln(10^20).")
print(f"For this value of t, we have e^t = {et} and e^(-t) = {e_neg_t}.")
print("\nPlugging these numbers into the solution equation gives:")
print(f"x(ln(10^20)) = (2 * (2 - {e_neg_t})) / ({et} + {e_neg_t})")
print(f"\nThe calculated value is:")
print(result)