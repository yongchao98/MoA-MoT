import sympy

# Define the symbolic variables used in the problem statement.
# y: vertical coordinate
# H: height of the microchannel
# k: Debye-Huckel parameter
# beta: slip length
# z1: zeta potential at the bottom surface (y=0)
y, H, k, beta, z1 = sympy.symbols('y H k beta z_1')

# The final expression for the potential psi(y) is derived from the
# linearized Poisson-Boltzmann equation subject to the given boundary conditions:
# BC1: psi(y=0) = z_1*(1 + beta*k)
# BC2: psi(y=H) = 0

# The problem asks to output each component of the final equation. We will construct
# the equation from its constituent parts.

# Part 1: The expression for the slip-dependent zeta potential at y=0.
zeta_potential_part = z1 * (1 + beta * k)

# Part 2: The numerator of the spatial distribution term, which is a hyperbolic sine function.
numerator_part = f"sinh(k*(H - y))"

# Part 3: The denominator of the spatial distribution term, which normalizes the expression.
denominator_part = f"sinh(k*H)"

# Now, we combine these parts into the final equation for the EDL potential psi(y).
# The following print statement displays the full, derived expression.
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
print(f"psi(y) = ({zeta_potential_part}) * {numerator_part} / {denominator_part}")