import math

# Let m and n be the dimensions of the source and target simplexes, respectively.
m = 23
n = 37

# The number of internal adjunctions from [m] to [n] in the simplex category
# is given by the binomial coefficient C(n + m, m).
# In our case, this is C(37 + 23, 23) = C(60, 23).
total = n + m
k = m
result = math.comb(total, k)

# Print the final equation and its result. The numbers in the final equation
# are 60 and 23.
print(f"C({total}, {k}) = {result}")