import sympy

# --- Plan Summary ---
# We are finding the Wilson-Fisher fixed point for phi^4 theory.
# This requires solving the equation for the one-loop beta function, beta(u*) = 0.
# The beta function is: beta(u) = -epsilon*u + (3*u^2)/(16*pi^2).
# We solve for the non-trivial fixed point u*.

# Define symbolic variables for the equation
u_star, epsilon = sympy.symbols('u^* epsilon')

# Set up the equation beta(u*) = 0, assuming u* is not zero
# -epsilon + 3*u* / (16*pi**2) = 0
equation = sympy.Eq(-epsilon + 3 * u_star / (16 * sympy.pi**2), 0)

# Solve for u*
fixed_point_solution = sympy.solve(equation, u_star)[0]

# --- Final Expression ---
# The prompt asks to output each number in the final equation.
# The expression is u* = (16 * pi^2 / 3) * epsilon.
# We will print this expression in a formatted string.
numerator_val = 16
exponent_val = 2
denominator_val = 3

print("The leading order expression for the fixed point coupling u^* in the phi^4 theory is:")
print(f"u^* = ({numerator_val} * \u03c0^{exponent_val} / {denominator_val}) * \u03f5")
