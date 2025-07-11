import math

# Number of terms in the sequence
N = 100000

# The necessary condition is M >= log10(N)
# We need to find the smallest positive integer M that satisfies this.

M_min = math.log10(N)

# The condition is M >= M_min.
# Since M must be an integer, the smallest integer M is ceil(M_min).
# However, the problem asks for the smallest *positive* integer.
# If M_min is, for example, 0.3, the smallest integer is 1. If M_min is 5.0, smallest integer is 5.

smallest_integer_M = math.ceil(M_min) if M_min > 0 else 1
if M_min == int(M_min) and M_min > 0: # handle exact integer cases
    smallest_integer_M = int(M_min)

print(f"Let N = {N}")
print("The condition for M to satisfy is M >= log10(N).")
print(f"M >= log10({N})")
print(f"M >= {M_min}")
print(f"The smallest positive integer M that satisfies this condition is {smallest_integer_M}.")