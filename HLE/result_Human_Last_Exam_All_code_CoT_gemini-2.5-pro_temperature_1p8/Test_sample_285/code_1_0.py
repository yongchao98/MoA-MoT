from fractions import Fraction

# The problem is to find the largest p such that the function I is not in L^p(R^9).
# This critical exponent, p, is determined by the dimension of the parameter space, D,
# and the slowest decay rate of the integral, delta_min.

# Step 1: Define the dimension of the parameter space.
# The integral I is a function of nine parameters (a_1, ..., a_9).
# The function I(a) is defined on R^9, so the dimension D is 9.
D = 9

# Step 2: Determine the slowest decay rate (delta_min) of the integral I(a).
# The decay rate of the oscillatory integral I(a) as |a| -> infinity is given by
# |I(a)| ~ |a|^(-delta). The rate delta depends on the geometry of the phase polynomial.
# To find the slowest decay (smallest delta), we must find the choice of coefficients 'a'
# that makes the phase maximally degenerate.
# For the given polynomial phase, this occurs when it reduces to a perfect cube of a
# linear form, for instance, A*x^3. For such a phase, the integral containing it decays
# with a rate of 1/3 according to van der Corput's Lemma.
# Thus, the minimum possible decay rate is delta_min = 1/3.
delta_min = Fraction(1, 3)

# Step 3: Calculate the critical exponent p.
# The condition for the L^p integral of I(a) to diverge is p * delta_min <= D.
# We are looking for the largest p for which I is not in L^p, which corresponds to
# the boundary case where the integral just diverges (p * delta_min = D).
p = D / delta_min

# Step 4: Output the reasoning and the final answer.
print("The largest p such that I is not in L^p(R^9) is given by the formula p = D / delta_min.")
print(f"The dimension of the parameter space is D = {D}.")
print(f"The slowest decay rate of the integral is delta_min = {delta_min.numerator}/{delta_min.denominator}.")
print(f"The final equation is p = {D} / ({delta_min.numerator}/{delta_min.denominator}).")
# Using the values D=9 and delta_min=1/3:
final_p = D / delta_min
print(f"p = {D} * {delta_min.denominator} / {delta_min.numerator} = {final_p:.0f}")
