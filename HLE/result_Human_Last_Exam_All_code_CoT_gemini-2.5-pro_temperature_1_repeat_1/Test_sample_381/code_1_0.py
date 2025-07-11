import math

# The problem is to find an upper-bound for ||B Q_{0, M}||_inf.
# The final form of the bound is C * sqrt(N), where N is the number of nodes
# and C is a constant independent of N.

# The question asks for the upper-bound expressed as a factor of sqrt(N).
# This factor is the constant C.
# The value of C depends on other parameters (epsilon, c, beta) which are not
# given numerically. Therefore, we represent it symbolically as "C".

# The prompt requires outputting each number in the final equation.
# The equation for the factor is simply "factor = C". We can write this as
# "factor = 1 * C" to include a number.
coefficient = 1
symbol = "C"

print("The upper bound for ||B Q_{0, M}||_inf is of the form C * sqrt(N).")
print(f"Here, C is a constant that is independent of N.")
print("The question asks for this upper bound expressed as a factor of sqrt(N).")
print("This factor is the constant C.")
print("The final equation for the factor is:")
print(f"{coefficient} * {symbol}")