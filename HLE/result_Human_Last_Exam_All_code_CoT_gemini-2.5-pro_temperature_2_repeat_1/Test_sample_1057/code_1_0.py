import sympy

# Define the symbols
pi, epsilon_0, a, V = sympy.symbols('pi epsilon_0 a V')

# Step 1: Define initial properties
# Initial radius = a
# Initial potential = V

# Step 2: Calculate initial capacitance and charge
# Capacitance of a sphere C = 4 * pi * epsilon_0 * r
C_initial = 4 * pi * epsilon_0 * a
# Initial charge Q = C * V
Q_initial = C_initial * V

# Step 3: Calculate the initial electrostatic energy
# Electrostatic energy U = 1/2 * Q * V
U_initial = sympy.Rational(1, 2) * Q_initial * V

# Step 4: Final energy is zero as the charge leaks completely
U_final = 0

# Step 5: The total dissipated Joule heat is the change in electrostatic energy
H = U_initial - U_final

# Step 6: Print the result.
# The problem asks to output the numbers in the final equation.
# The derived expression for H is 2 * pi * epsilon_0 * a * V**2.
# The numbers are 2 and the implicit exponent 2.

# For a clearer output, let's substitute the numeric value for the fraction.
H_simplified = sympy.simplify(H)
# Extract coefficients and variables to print them cleanly
coeff, terms = H_simplified.as_coeff_terms()

print("The formula for the total Joule heat (H) dissipated is derived from the initial stored electrostatic energy.")
print(f"The initial energy is U_initial = (1/2) * Q_initial * V.")
print(f"With Q_initial = (4 * pi * epsilon_0 * a) * V, this becomes:")
print(f"U_initial = (1/2) * (4 * pi * epsilon_0 * a * V) * V = 2 * pi * epsilon_0 * a * V^2.")
print("\nSince all of this initial energy is ultimately dissipated as heat, the total Joule heat is:")
# Manually constructing the print statement to show the numbers
print("H = 2 * pi * epsilon_0 * a * V**2")