# The dimension of the space is n.
# The moment curve is in R^3, so n=3.
n = 3

# The critical exponent p_c is given by the formula n*(n+1)/2.
# For p > p_c, there exist non-zero L^p functions whose Fourier transform
# is supported on the moment curve.
# For p <= p_c, any L^p function with Fourier support on the curve must be zero.
# The question asks for the largest possible value of p for which this holds, which is p_c.
p_c = n * (n + 1) / 2

# We print the final equation and the result.
print(f"The largest possible value of p is calculated using the formula p = n * (n + 1) / 2.")
print(f"For n = {n}, the calculation is:")
# The final output needs to show the numbers in the equation.
print(f"p = {n} * ({n} + 1) / 2 = {int(p_c)}")