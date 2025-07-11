import math

# Plan:
# The constant `b` is determined by combining two key theorems about random walks.
# 1. The relationship between the expected cover-and-return time C(n) and the
#    expected cover time E[T_cov]: C(n) = 2 * E[T_cov].
# 2. The asymptotic formula for the expected cover time on a random tree:
#    E[T_cov] ~ sqrt(pi / 2) * n^(3/2).
#
# Combining them gives:
# C(n) ~ 2 * sqrt(pi / 2) * n^(3/2)
# C(n) ~ sqrt(4 * pi / 2) * n^(3/2)
# C(n) ~ sqrt(2 * pi) * n^(3/2)
#
# The problem asks for the constant `b` in C(n) ~ b * n^(3/2).
# Thus, b = sqrt(2 * pi).

# The final equation for b is b = sqrt(2 * pi).
# Here are the numbers used in that equation.
val_2 = 2
val_pi = math.pi

b = math.sqrt(val_2 * val_pi)

print("The formula for the constant b is derived as follows:")
print("C(n) = E[cover-and-return time] = 2 * E[cover time]")
print("E[cover time] is asymptotic to sqrt(pi / 2) * n^(3/2)")
print("So, b = 2 * sqrt(pi / 2) = sqrt(4 * pi / 2) = sqrt(2 * pi)")
print("\n--- Final Equation ---")
print("b = sqrt(val_2 * val_pi)")
print(f"Value of 'val_2': {val_2}")
print(f"Value of 'val_pi': {val_pi}")
print("\n--- Final Answer ---")
print("The exact value of b is sqrt(2 * pi).")
print("The numerical value of b is:")
print(b)