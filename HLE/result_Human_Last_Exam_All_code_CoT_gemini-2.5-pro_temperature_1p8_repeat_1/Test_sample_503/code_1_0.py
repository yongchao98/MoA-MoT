import math

# The problem is to find the number of internal adjunctions in the simplex
# category from [n] to [m]. This is given by the binomial coefficient C(n+m, n).
n = 23
m = 37

# The values for the binomial coefficient C(a, b)
a = n + m
b = n

# Calculate the result
result = math.comb(a, b)

# The final equation is C(60, 23) = result.
# The instruction is to output each number in this final equation.
# The following print statement displays the equation and all its numbers.
print(f"C({a}, {b}) = {result}")
