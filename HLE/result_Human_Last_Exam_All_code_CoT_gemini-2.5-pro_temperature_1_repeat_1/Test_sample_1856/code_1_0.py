import math

# Step 1: Define the given velocity of particle A.
beta_A = 0.95

# Step 2: Calculate the Lorentz factor, gamma_A.
gamma_A = 1 / math.sqrt(1 - beta_A**2)

# Step 3: Calculate the tangent of the angle theta_C in the lab frame.
# The formula is tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A)).
tan_theta_C = 1 / (gamma_A * (1 + math.sqrt(2) * beta_A))

# Step 4: Calculate the angle theta_C in degrees from its tangent.
theta_C_deg = math.degrees(math.atan(tan_theta_C))

# Step 5: Print the final equation with the numbers and the result.
# The instruction is to output each number in the final equation.
print("The final equation for the angle theta_C is:")
print(f"theta_C = arctan(1 / (gamma_A * (1 + sqrt(2) * beta_A)))")
print("\nSubstituting the numerical values:")
# We show the final calculation in one line with intermediate values for clarity.
final_calculation = f"theta_C = arctan(1 / ({gamma_A:.3f} * (1 + {math.sqrt(2):.3f} * {beta_A}))) = {theta_C_deg:.3f} degrees"
print(final_calculation)