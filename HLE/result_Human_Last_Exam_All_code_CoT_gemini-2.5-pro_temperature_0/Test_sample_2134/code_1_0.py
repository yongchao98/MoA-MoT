import math

# Step 1: Define the time t
t = 2 * math.sqrt(3)

# Step 2: Define the components of the initial position x(0)
# x(0) = 3 + (6*(3-sqrt(3)))**(1/3) + (6*(3+sqrt(3)))**(1/3)
# Let's define the terms inside the cube roots
sqrt3 = math.sqrt(3)
term_A_base = 18 - 6 * sqrt3
term_B_base = 6 * (3 + sqrt3) # This is the same as 18 + 6*sqrt(3)

# Step 3: The trajectory equation is x(t) = x(0) - t^2 / 4
# Calculate the displacement term t^2 / 4
displacement = t**2 / 4

# Step 4: Calculate the final position x(t)
# x(t) = (3 + term_A_base**(1/3) + term_B_base**(1/3)) - displacement
# Since displacement is 3, this simplifies to:
# x(t) = term_A_base**(1/3) + term_B_base**(1/3)
term_A = term_A_base**(1/3)
term_B = term_B_base**(1/3)
final_position = term_A + term_B

# Step 5: Print the equation and the result, showing each number.
print(f"The final position is calculated using the formula: x(t) = x(0) - t^2 / 4")
print(f"Given x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)")
print(f"Given t = 2*sqrt(3)")
print(f"The displacement is t^2 / 4 = ({t:.4f})^2 / 4 = {t**2:.4f} / 4 = {displacement:.4f}")
print(f"So, x({t:.4f}) = x(0) - {displacement:.4f}")
print(f"x({t:.4f}) = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)")
print(f"Let's evaluate the terms:")
print(f"First term: ({term_A_base:.4f})^(1/3) = {term_A:.4f}")
print(f"Second term: ({term_B_base:.4f})^(1/3) = {term_B:.4f}")
print(f"Final position = {term_A:.4f} + {term_B:.4f} = {final_position:.4f}")
print(f"The exact value is the real root of the cubic equation y^3 - 18*y - 36 = 0.")
print(f"The final numerical value is {final_position}")
