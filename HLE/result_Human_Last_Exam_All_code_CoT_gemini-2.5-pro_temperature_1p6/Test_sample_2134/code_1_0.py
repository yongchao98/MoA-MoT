import numpy as np

# Step 1: Define and calculate the initial position x(0)
# x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
sqrt3 = np.sqrt(3)

# Define the terms in the x(0) expression
term1_x0 = 3
term2_x0 = (18 - 6 * sqrt3)**(1/3)
term3_x0 = (18 + 6 * sqrt3)**(1/3)

# Calculate x(0)
x0 = term1_x0 + term2_x0 + term3_x0

# Step 2: Define the time t
t = 2 * sqrt3

# Step 3: Calculate the displacement, which is t^2 / 4
displacement = t**2 / 4

# Step 4: Calculate the final position x(t) using the classical trajectory formula
# x(t) = x(0) + t^2 / 4
xt = x0 + displacement

# Step 5: Print the detailed calculation
print("Calculation Steps:")
print(f"The initial position is x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)")
print(f"Term 2 of x(0) = (18 - 6 * {sqrt3:.4f})^(1/3) = ({18 - 6 * sqrt3:.4f})^(1/3) = {term2_x0:.4f}")
print(f"Term 3 of x(0) = (18 + 6 * {sqrt3:.4f})^(1/3) = ({18 + 6 * sqrt3:.4f})^(1/3) = {term3_x0:.4f}")
print(f"So, x(0) = {term1_x0} + {term2_x0:.4f} + {term3_x0:.4f} = {x0:.4f}\n")

print(f"The particle evolves for a time t = 2 * sqrt(3) = {t:.4f}")
print(f"The displacement is given by t^2 / 4 = ({t:.4f})^2 / 4 = {displacement:.4f}\n")

print("The final position is x(t) = x(0) + t^2 / 4")
print(f"x({t:.4f}) = {x0:.4f} + {displacement:.4f}")
print(f"x({t:.4f}) = {xt:.4f}")

# Final Answer
print("\nFinal equation with evaluated numbers:")
print(f"x(2*sqrt(3)) = {x0:.4f} + ({t:.4f})^2 / 4 = {xt:.4f}")
print(f"{xt:.4f}")
<<<11.0171>>>