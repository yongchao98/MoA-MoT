import math

# Step 1: Define initial parameters
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Step 2: Calculate the constant k from R = k * E^(3/2)
# k = R0 / E0^(3/2)
k = R0 / (E0 ** 1.5)
print(f"The Geiger's rule constant for air (k) is calculated as: {R0} / ({E0}^1.5) = {k:.4f} cm/MeV^1.5\n")

# Step 3: Calculate the energy E(x) at distance x
# The residual range at x=4 cm is R_residual = R0 - x
R_residual = R0 - x
# From R_residual = k * E(x)^(3/2), we can find E(x)
# E(x) = (R_residual / k)^(2/3)
E_at_x = (R_residual / k) ** (2.0 / 3.0)
print(f"The residual range at {x} cm is: {R0} - {x} = {R_residual:.1f} cm")
print(f"The energy of the α-particle at {x} cm from the source is E({x}) = ({R_residual:.1f} / {k:.4f})^(2/3) = {E_at_x:.4f} MeV\n")

# Step 4: Calculate the energy loss per centimeter (-dE/dx) at x
# We have R = k * E^(3/2).
# Differentiating with respect to E: dR/dE = (3/2) * k * E^(1/2)
# Since dx = -dR, the energy loss is -dE/dx = dE/dR = 1 / (dR/dE)
# -dE/dx = 1 / ((3/2) * k * sqrt(E))
energy_loss = 1 / (1.5 * k * math.sqrt(E_at_x))

print(f"The energy loss per centimeter (-dE/dx) is calculated using the formula: 1 / (1.5 * k * sqrt(E({x})))")
print(f"Plugging in the values:")
print(f"-dE/dx = 1 / (1.5 * {k:.4f} * sqrt({E_at_x:.4f}))")
print(f"-dE/dx = {energy_loss:.4f} MeV/cm")
<<<0.8540>>>