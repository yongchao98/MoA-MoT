# Initial parameters given in the problem
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Based on Geiger's rule, the energy loss per centimeter (-dE/dx) is calculated.
# The formula is: (2/3) * E0 * R0^(-2/3) * (R0 - x)^(-1/3)
energy_loss_per_cm = (2/3) * E0 * (R0**(-2/3)) * ((R0 - x)**(-1/3))

# Print the final equation with the numerical values and the result.
# This fulfills the request to show each number in the final equation.
print(f"Energy loss calculation (-dE/dx):")
print(f"(2/3) * {E0} * ({R0})**(-2/3) * ({R0} - {x})**(-1/3) = {energy_loss_per_cm:.3f} MeV/cm")