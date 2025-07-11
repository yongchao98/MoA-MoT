import sympy

# Define symbolic variables for resistance and capacitance
r = sympy.Symbol('r')
C = sympy.Symbol('C')

# The equivalent Thevenin resistance R_th at node a0 for the
# infinite ladder network is found by solving the circuit's difference equations.
# The derivation shows that R_th = r * (1 + sqrt(3)).
# For completeness, we can derive this result symbolically.
k = sympy.Symbol('k')
# The resistance calculation depends on the root of the characteristic equation k^2 - 4*k + 1 = 0
# that has a magnitude less than 1.
roots = sympy.solve(k**2 - 4*k + 1, k)
k1 = roots[0]  # This is 2 - sqrt(3)

# The Thevenin resistance R_th is derived to be 2*r / (1 - k1).
R_th = sympy.simplify(2*r / (1 - k1))   # This evaluates to r*(1 + sqrt(3))

# The cutoff frequency omega_c is given by the formula omega_c = 1 / (R_th * C).
omega_c = 1 / (R_th * C)

# To fulfill the request to output each number in the final equation,
# we can deconstruct the symbolic expression.
# The final expression is omega_c = 1 / (r*C*(1 + sqrt(3)))
numerator = 1
# The coefficient of r*C in the denominator is (1 + sqrt(3))
denominator_coefficient = R_th / r
# This coefficient is an addition of two terms.
term_1 = denominator_coefficient.args[0]
sqrt_term = denominator_coefficient.args[1]
number_under_sqrt = sqrt_term.args[0]

# Print the final equation for the cutoff frequency, showing each numerical component.
print("The cutoff frequency at node a0 is:")
print(f"omega_c = {numerator} / (({term_1} + sqrt({number_under_sqrt})) * r * C)")
