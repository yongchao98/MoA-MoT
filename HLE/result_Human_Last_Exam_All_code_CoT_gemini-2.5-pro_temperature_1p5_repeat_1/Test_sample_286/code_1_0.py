import math

# The condition derived from the analysis is M >= log10(2).
# We need to find the smallest positive integer M that satisfies this.

log10_2 = math.log10(2)
# The condition is M >= log10_2.
# Since M must be a positive integer, we need to find the smallest integer greater than or equal to log10_2.
# This is equivalent to taking the ceiling of log10_2.
M = math.ceil(log10_2)

# However, the problem specifies M must be a positive integer.
# math.ceil(0.301...) = 1, which is a positive integer.
# So, the smallest positive integer M is 1.

print(f"The condition is M >= log10(2)")
print(f"log10(2) is approximately {log10_2}")
print(f"The smallest positive integer M that satisfies this condition is {M}")
