import math

# Step 1: Define the solution from our derivation.
# We found that c^2 = 5/2 at the point of emergence of the giant component.
c_squared_num = 5
c_squared_den = 2
c_squared = c_squared_num / c_squared_den

# The exact value of c is the square root of this value.
c = math.sqrt(c_squared)
c_fourth = c_squared**2

# Step 2: Define the terms of the criticality equation E[d^2]/E[d] = 2.
# The equation is: (c^2/3 + 2*c^4/15) / (c^2/3) = 2
numerator_term1 = c_squared / 3
numerator_term2 = 2 * c_fourth / 15
denominator_term = c_squared / 3

# Step 3: Print the verification of the equation.
# This shows how each number in the equation is calculated from the value of c.
print("The criticality condition for the emergence of the giant component is E[d^2] / E[d] = 2.")
print("Based on the model, we derived E[d] = c^2/3 and E[d^2] = c^2/3 + 2*c^4/15.")
print(f"Solving this gives c^2 = {c_squared_num}/{c_squared_den} = {c_squared}.\n")
print("We verify this by plugging c back into the equation.")

print("The numbers in the final equation are:")
print(f"Numerator term 1 (c^2/3): {numerator_term1:.4f}")
print(f"Numerator term 2 (2*c^4/15): {numerator_term2:.4f}")
print(f"Denominator term (c^2/3): {denominator_term:.4f}\n")

numerator_val = numerator_term1 + numerator_term2
denominator_val = denominator_term
kappa = numerator_val / denominator_val

print("The final equation is:")
print(f"({numerator_term1:.4f} + {numerator_term2:.4f}) / {denominator_term:.4f} = {kappa}")
print(f"This simplifies to {numerator_val:.4f} / {denominator_val:.4f} = {kappa}, which confirms the condition.\n")

# Step 4: Print the final answer.
print(f"The exact value of the constant c is the square root of 5/2.")
print(f"c = {c}")