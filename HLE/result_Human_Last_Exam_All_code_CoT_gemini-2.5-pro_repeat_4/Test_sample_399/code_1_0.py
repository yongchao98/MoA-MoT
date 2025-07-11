from fractions import Fraction

# We are solving for x(1) from the particular solution:
# -1/x(t) = t^2/2 - t + 1/8

# Set the value of t for which we want to find x
t = 1

# The constant C was found to be 1/8 from the initial condition x(0) = -8
C = 1/8

# The final equation involves substituting t=1. We will print the numbers
# involved in this equation.
# Equation: -1/x(1) = (1**2)/2 - 1 + 1/8
term1_num = 1**2
term1_den = 2
term2 = -1
term3_num = 1
term3_den = 8

print("The final equation to solve is derived by substituting t=1 into the particular solution:")
print(f"-1/x(1) = ({term1_num**2})/{term1_den} + ({term2}) + {term3_num}/{term3_den}")

# Calculate the value of the right-hand side (RHS) of the equation
rhs = (t**2 / 2) - t + C

# Calculate the final value of x(1)
x_1 = -1 / rhs

# Print the intermediate and final results
print(f"\nCalculating the right-hand side:")
print(f"-1/x(1) = {rhs}")

# For a more readable exact answer, we convert the result to a fraction
final_answer_fraction = Fraction(x_1).limit_denominator()

print(f"\nSolving for x(1):")
print(f"x(1) = {final_answer_fraction}")