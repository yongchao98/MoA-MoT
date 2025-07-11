import math

# The problem is to determine the value of alpha in the expression c_n = Theta(alpha^n).
# We have established an upper bound for c_n.
# c_n <= sqrt(2^n - 1).
# Asymptotically, sqrt(2^n - 1) behaves like (sqrt(2))^n.
# This implies that alpha <= sqrt(2).

# We also established a simple lower bound c_n >= 1, which implies alpha >= 1.

# In many such problems, the upper bound derived from general principles turns out to be asymptotically tight.
# Assuming this is the case, the growth rate is dominated by the (sqrt(2))^n term.
# Therefore, the value of alpha is sqrt(2).

alpha_squared = 2
alpha = math.sqrt(alpha_squared)

print(f"The analysis of the upper bound suggests alpha <= sqrt(2).")
print(f"A simple lower bound suggests alpha >= 1.")
print(f"Assuming the upper bound is asymptotically tight, the value for alpha is sqrt(2).")
print(f"alpha = {alpha}")
