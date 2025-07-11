import sympy

# Define the variable for the squared mass given in the problem
m_squared = sympy.Symbol('m^2')

# From the analysis of the Lagrangian for the scalar (trace) mode h,
# we have the general form L = (A/2)*(dh)^2 - (B/2)*h^2.
# The squared mass is then M^2 = B/A.

# Coefficient of the kinetic term for h relative to the canonically normalized spin-2 part.
# This comes from the structure of the linearized Einstein-Hilbert action.
# The Lagrangian for the trace mode has a relative factor of -3/4 compared to the traceless mode.
# L_kin = 1/2 (d h_bar)^2 - 3/8 (d h)^2.
# So, the kinetic term is L_kin_scalar = -3/8 * (dh)^2 = (-3/4)/2 * (dh)^2.
A = sympy.Rational(-3, 4)

# The potential term for the scalar part is V_scalar = (m^2/2) * (1/4 * h^2) = (m^2/8) * h^2.
# This corresponds to L_mass_scalar = -V_scalar = -(m^2/8) * h^2 = -(m^2/4)/2 * h^2
B = m_squared / 4

# Calculate the squared mass of the scalar mode
M_squared = B / A

print("The Lagrangian for the scalar degree of freedom (h) is of the form L = A/2 * (∂h)^2 - B/2 * h^2")
print(f"The kinetic term coefficient A is: {A}")
print(f"The mass term coefficient B is: {B}")
print("\nThe squared mass of the sixth degree of freedom is B/A:")
# Use pretty print for the final equation
equation = sympy.Eq(sympy.Symbol('M^2'), M_squared, evaluate=False)
sympy.pprint(equation, use_unicode=False)
print(f"\nFinal calculated value: {M_squared}")