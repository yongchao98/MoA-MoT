import math

# This script calculates the minimum probability q (q_0) for Alice to win with at least 95% probability.

# Alice's winning probability, p_A
p_A = 0.95

# The relationship between p_A and q_0 is given by the equation:
# p_A = 1 - (1 - q_0**4 * p_A**3)**3
# We need to solve for q_0 when p_A = 0.95.

# Rearranging the equation to solve for q_0:
# q_0 = ((1 - (1 - p_A)**(1/3)) / p_A**3)**(1/4)

# Let's break down the calculation and display the numbers involved in the equation.
# Numbers from the equation:
val_p_A = p_A
val_one_minus_p_A = 1 - p_A
exp_one_third = 1.0 / 3.0
exp_three = 3.0
exp_one_fourth = 1.0 / 4.0

print(f"We are solving the equation: {val_p_A} = 1 - (1 - q_0**4 * {val_p_A}**{exp_three})**{exp_three}")
print("\nThe rearranged equation for q_0 is:")
print(f"q_0 = ((1 - (1 - {val_p_A})**({exp_one_third:.4f})) / {val_p_A}**{exp_three})**({exp_one_fourth:.4f})")
print(f"q_0 = ((1 - {val_one_minus_p_A:.2f}**({exp_one_third:.4f})) / {val_p_A**exp_three:.8f})**({exp_one_fourth:.4f})")

# Perform the calculation
numerator = 1 - (1 - p_A)**(1/3)
denominator = p_A**3
q_0 = (numerator / denominator)**(1/4)

print(f"\nCalculated value of q_0: {q_0:.8f}")

# The problem asks for the floor of 100 * q_0
result = math.floor(100 * q_0)

print(f"\nThe value of 100 * q_0 is: {100 * q_0:.8f}")
print(f"The final answer, floor(100 * q_0), is: {result}")
