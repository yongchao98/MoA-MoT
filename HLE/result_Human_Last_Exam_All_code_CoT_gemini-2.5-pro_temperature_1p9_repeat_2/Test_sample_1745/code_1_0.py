import sympy

# This script finds the expression for the Electrical Double-Layer (EDL) potential distribution
# in a parallel-plate microchannel based on the linearized Poisson-Boltzmann equation.

# Define the mathematical symbols required for the equation.
# y: the spatial coordinate perpendicular to the channel walls
# H: the height of the microchannel
# k: the Debye-Hückel parameter
# z_1: the zeta potential at the bottom surface (y=0)
# psi: the function representing the EDL potential
y, H, k = sympy.symbols('y H k')
z_1 = sympy.Symbol('z_1')
psi = sympy.Function('psi')

# The expression for the EDL potential distribution is the solution to the
# linearized Poisson-Boltzmann equation (d^2(psi)/dy^2 = k^2 * psi)
# with the boundary conditions: psi(0) = z_1 and psi(H) = 0.
final_expression = z_1 * sympy.sinh(k * (H - y)) / sympy.sinh(k * H)

# Create a sympy Equality object to represent the final equation.
final_equation = sympy.Eq(psi(y), final_expression)

# Print the final equation. The output string shows all the components
# of the equation, including the variables and mathematical functions.
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
print(str(final_equation))
