# My analysis concluded that for both regimes, the optimal query complexity for sorting is Theta(N * log N).
# To convert this into the requested (a,b,c) notation, we solve the equation:
# (N * log N)^2 = N^a * (log N)^b * (log log N)^c
#
# which can be written as:
# N^2 * (log N)^2 * (log log N)^0 = N^a * (log N)^b * (log log N)^c
#
# By comparing the exponents of the corresponding terms on both sides of the equation,
# we determine the integer values for a, b, and c.

# For the first regime, the parameters are:
a1 = 2  # from the exponent of N
b1 = 2  # from the exponent of log N
c1 = 0  # from the exponent of log log N

# For the second regime, the analysis also yields a complexity of Theta(N * log N).
# Therefore, the parameters are the same.
a2 = 2
b2 = 2
c2 = 0

# The final answer is presented in the format (a,b,c),(a,b,c)
print(f"({a1},{b1},{c1}),({a2},{b2},{c2})")