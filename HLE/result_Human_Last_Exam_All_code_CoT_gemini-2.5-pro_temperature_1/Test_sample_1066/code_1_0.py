import math

# The derived lower bound for E[S] is n*alpha + n*C
# where C is the constant 2 - pi^2/6.
# This code calculates and prints the numbers in the final equation's constant term.

pi_squared = math.pi**2
constant_term_per_token = 2 - pi_squared / 6
coefficient_of_n_alpha = 1

print("The lower bound for E[S] is given by the equation: E[S] >= n*alpha + n*C")
print(f"The coefficient of n*alpha is: {coefficient_of_n_alpha}")
print(f"The constant C is 2 - pi^2/6, which evaluates to:")
print(constant_term_per_token)

# The final equation is E[S] >= 1 * n*alpha + n*(2 - pi^2/6)
# The numbers in the equation are 1, 2, 6, and pi.
# The prompt asks to output each number in the final equation.
# The numbers defining the constant C are 2 and 6.
print("The integer numbers defining the constant term are 2 and 6.")
