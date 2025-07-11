# c represents the cardinality of the continuum.
# The expression for the largest possible weight is 2^(2^c).
# This code prints the components of this expression.

power1_base = 2
power1_exp_base = 2
power1_exp_exp = "c"

print(f"The largest possible weight is {power1_base}^({power1_exp_base}^{power1_exp_exp})")