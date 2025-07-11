import math

# Define the initial parameters from the problem
initial_energy_E0 = 8.5  # in MeV
total_range_R0 = 8.3     # in cm
distance_x = 4.0         # in cm

# Calculate the numerator of the formula
numerator = 2 * initial_energy_E0

# Calculate the denominator of the formula
remaining_range = total_range_R0 - distance_x
denominator = 3 * (total_range_R0**(2/3)) * (remaining_range**(1/3))

# Calculate the energy loss per centimeter (-dE/dx)
energy_loss_per_cm = numerator / denominator

# Print the explanation and the result
print("To solve this problem, we use a formula derived from Geiger's rule for the range of alpha particles.")
print("The formula for energy loss per centimeter (-dE/dx) at a distance x is:")
print("-dE/dx = (2 * E₀) / (3 * R₀^(2/3) * (R₀ - x)^(1/3))\n")

print("Plugging in the given values:")
print(f"Initial Energy (E₀) = {initial_energy_E0} MeV")
print(f"Total Range (R₀) = {total_range_R0} cm")
print(f"Distance (x) = {distance_x} cm\n")

print("The calculation is:")
print(f"-dE/dx = (2 * {initial_energy_E0}) / (3 * ({total_range_R0})^(2/3) * ({total_range_R0} - {distance_x})^(1/3))")
print(f"-dE/dx = {numerator} / (3 * {total_range_R0**(2/3):.4f} * {remaining_range**(1/3):.4f})")
print(f"-dE/dx = {numerator} / {denominator:.4f}\n")

print(f"The calculated energy loss per centimetre at {distance_x} cm is: {energy_loss_per_cm:.3f} MeV/cm")
<<<0.849>>>