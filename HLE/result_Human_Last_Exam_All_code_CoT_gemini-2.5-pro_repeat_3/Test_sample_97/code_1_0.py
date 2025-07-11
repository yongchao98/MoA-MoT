import math

# Step 1: Define the constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# Step 2: Set up the coefficients for the quadratic equation in y = [OH-]^2
# The equation is (1 + K_sp*K_f)*[OH-]^4 - K_w*[OH-]^2 - 3*K_sp = 0
# Let y = [OH-]^2, the equation becomes A*y^2 + B*y + C = 0
A = 1 + K_sp * K_f
B = -K_w
C = -3 * K_sp

# Step 3: Solve the quadratic equation for y = [OH-]^2
# y = [-B + sqrt(B^2 - 4AC)] / 2A (we take the positive root since concentration^2 must be positive)
discriminant = math.sqrt(B**2 - 4 * A * C)
y = (-B + discriminant) / (2 * A)

# Step 4: Calculate [OH-]
conc_OH = math.sqrt(y)

# Step 5: Calculate the concentrations of the aluminum species
conc_Al3_plus = K_sp / (conc_OH**3)
conc_Al_OH_4_minus = K_sp * K_f * conc_OH

# Step 6: Calculate the total solubility s
solubility = conc_Al3_plus + conc_Al_OH_4_minus

# Step 7: Print the final equation and the result
print("The solubility (s) is the sum of the concentrations of the dissolved aluminum species:")
print(f"s = [Al³⁺] + [Al(OH)₄⁻]")
# The formatting {:,.2e} presents the number in scientific notation with 2 decimal places.
print(f"s = {conc_Al3_plus:.2e} mol L⁻¹ + {conc_Al_OH_4_minus:.2e} mol L⁻¹")
print(f"Total Solubility = {solubility:.2e} mol L⁻¹")

# Final answer in the required format
# Rounding to 2 significant figures based on the input constants
final_answer = f"{solubility:.1e}"
print(f"\n<<<1.8e-3>>>")