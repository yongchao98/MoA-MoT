import sympy

# This script calculates the leading order expression for the Wilson-Fisher fixed point
# in phi^4 theory using the epsilon expansion.

# 1. Define symbolic variables for the coupling (u), the expansion parameter (epsilon), and pi.
u, epsilon, pi = sympy.symbols('u epsilon pi')

# 2. Define the one-loop beta function in d = 4 - epsilon dimensions.
# The term -epsilon*u arises from the canonical dimension of the coupling.
# The term +3*u**2/(16*pi**2) is the one-loop quantum correction.
beta_function = -epsilon * u + (3 * u**2) / (16 * pi**2)

print("The one-loop beta function is: beta(u) = -epsilon*u + 3*u**2 / (16*pi**2)")
print("-" * 60)

# 3. Find the fixed points u* by solving the equation beta(u*) = 0.
fixed_points = sympy.solve(beta_function, u)

# 4. The solver returns a list of solutions. We are interested in the non-trivial one.
# We filter out the trivial solution u=0.
u_star_nontrivial = None
for sol in fixed_points:
    if sol != 0:
        u_star_nontrivial = sol
        break

print("A fixed point u* is a value where beta(u*) = 0.")
print(f"Solving the equation '{beta_function} = 0' for u gives two solutions:")
print(f"1. Trivial (Gaussian) fixed point: u* = 0")
print(f"2. Non-trivial (Wilson-Fisher) fixed point: u* = {u_star_nontrivial}")
print("-" * 60)

# 5. Extract the components of the non-trivial solution to display the final equation clearly.
# The expression is of the form: (numerator/denominator) * pi**power * epsilon
coeff = u_star_nontrivial / (pi**2 * epsilon)
numerator, denominator = coeff.as_numer_denom()
power = sympy.degree(u_star_nontrivial, pi) # This will be 2

print("The leading order expression for the non-trivial fixed point coupling u* is:")
# Final print statement showing the equation with numbers.
print(f"u* = ({numerator} * pi^{power} / {denominator}) * epsilon")