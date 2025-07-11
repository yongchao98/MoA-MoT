import math

# The normalised density of the invariant measure for the map T(x) = 1/x^(1/2) mod 1
# is given by the function rho(x).

# The unnormalized density is proportional to 1 / (x + sqrt(x)).
# To normalize it, we compute the integral of this function from 0 to 1.
# The value of the integral is 2 * ln(2).
# So, the normalization constant is 1 / (2 * ln(2)).

c0 = 1
c1 = 2
c2 = 2
c3 = 1

# The final equation is rho(x) = c0 / (c1 * ln(c2) * (x + x^(c3/c2)))
# which simplifies to rho(x) = 1 / (2 * ln(2) * (x + sqrt(x)))

print("The normalised density of the invariant measure is rho(x) = c0 / (c1 * ln(c2) * (x + sqrt(x)))")
print(f"where c0 = {c0}, c1 = {c1}, c2 = {c2}.")
print("The final equation is:")
print(f"rho(x) = {c0} / ({c1} * ln({c2}) * (x + sqrt(x)))")

# We can also print the approximate value of the normalization constant
norm_constant = 1 / (2 * math.log(2))
print(f"\nThe numerical value of the normalization constant 1/(2*ln(2)) is approximately {norm_constant:.4f}.")