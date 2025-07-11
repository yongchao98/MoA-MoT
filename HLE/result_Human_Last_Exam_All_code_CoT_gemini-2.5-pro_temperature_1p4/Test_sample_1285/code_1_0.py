# (a) The general formula for the maximum number of roots is t(t-1)/2.

# (b) We calculate this value for t = 5.
t = 5

# The formula gives the degree of the polynomial part of the Wronskian.
# This corresponds to the maximum number of roots in the interval ]0, 1[.
max_roots = t * (t - 1) // 2

# We print the calculation as requested, showing each number in the equation.
print(f"For t = {t}, the maximum number of roots is:")
print(f"({t} * ({t} - 1)) / 2 = {max_roots}")
