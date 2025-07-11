# This script prints the analytical expression for the initial magnetization
# curve of a thin superconducting strip in a perpendicular magnetic field,
# based on the critical-state model.

# The parameters of the problem are treated as symbolic strings to
# construct the final formula.
a = 'a'   # Half-width of the rectangular cross-section
b = 'b'   # Half-thickness of the rectangular cross-section
Jc = 'Jc' # Constant critical current density
H = 'H'   # Applied magnetic field strength

# The derivation results in an equation for magnetization M as a function of H.
# The structure of the equation is: M(H) = -Prefactor * tanh^2(Argument)
# The final equation must contain all numerical constants from the derivation.

# 1. Construct the prefactor term: (Jc * a / 4)
# This part includes the number 4.
prefactor = f"({Jc} * {a} / 4)"

# 2. Construct the argument for the hyperbolic tangent function (tanh):
# The argument is (π * H) / (2 * b * Jc).
# This part includes the number 2 and the constant π.
tanh_argument_numerator = f"π * {H}"
tanh_argument_denominator = f"2 * {b} * {Jc}"
tanh_argument = f"({tanh_argument_numerator} / {tanh_argument_denominator})"

# 3. Assemble the complete equation string.
final_equation = f"M({H}) = -{prefactor} * tanh^2{tanh_argument}"

# 4. Print the final analytical expression.
print("The analytical expression for the initial magnetization curve M(H) is:")
print(final_equation)