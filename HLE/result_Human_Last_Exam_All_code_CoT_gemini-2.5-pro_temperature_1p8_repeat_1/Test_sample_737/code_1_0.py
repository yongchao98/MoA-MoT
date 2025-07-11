import math

# Step 1: Define Given Information
# Plate Bearing Test on Subgrade
P_plate = 30000  # N
d_plate = 305 / 1000  # m
a_plate = d_plate / 2  # m
delta_subgrade = 2460 / 1e6  # m

# Plate Bearing Test on Trial Pavement
h_trial = 300 / 1000  # m
delta_trial = 1080 / 1e6  # m

# Design Wheel Load
W_design_kg = 1.80 * 1000  # kg
P_design = W_design_kg * 9.81 # N
q_design = 600 * 1000  # N/m^2 (Pa)
delta_design_limit = 1.00 / 1000  # m

# Material Properties
mu = 0.5  # Poisson's ratio for all materials

print("--- Step 1: Calculate Subgrade Modulus (Es) ---")
# Using Boussinesq's equation for deflection under a flexible circular load on an elastic half-space.
# For mu=0.5, the formula simplifies to: Delta = 1.5 * P / (pi * a * Es)
# Rearranging for Es: Es = 1.5 * P / (pi * a * Delta)
Es = (1.5 * P_plate) / (math.pi * a_plate * delta_subgrade)
print(f"The calculated subgrade modulus Es (E2) is: {Es/1e6:.2f} MPa")
print(f"Es = (1.5 * {P_plate} N) / (pi * {a_plate:.4f} m * {delta_subgrade:.6f} m) = {Es:.0f} Pa\n")


print("--- Step 2: Analyze Trial Pavement Section ---")
# The deflection for a two-layer system is given by: Delta = (1.5 * q * a / Es) * F2
# First, calculate the contact pressure q_plate for the plate test.
q_plate = P_plate / (math.pi * a_plate**2)
# Then, calculate the deflection factor F2 for the trial section.
F2_trial = (delta_trial * Es) / (1.5 * q_plate * a_plate)
# The geometry of the trial section is defined by the h/a ratio.
h_over_a_trial = h_trial / a_plate
print(f"The geometric ratio for the trial section is h/a = {h_trial:.3f} m / {a_plate:.4f} m = {h_over_a_trial:.3f}")
print(f"The deflection factor F2 for the trial section is calculated as: {F2_trial:.4f}\n")


print("--- Step 3: Determine Design Requirements ---")
# Calculate the radius of the design wheel load.
# P_design = q_design * pi * a_design^2
a_design = math.sqrt(P_design / (math.pi * q_design))
print(f"The design load contact radius is a_design = sqrt({P_design:.0f} N / (pi * {q_design:.0f} Pa)) = {a_design:.4f} m")

# Calculate the required deflection factor F2 for the design case to meet the deflection limit.
F2_design = (delta_design_limit * Es) / (1.5 * q_design * a_design)
print(f"The required deflection factor F2 for the design is calculated as: {F2_design:.4f}\n")


print("--- Step 4: Determine Final Pavement Thickness ---")
print("Observation: The required design factor F2_design is approximately equal to the trial factor F2_trial.")
print(f"({F2_design:.4f} vs {F2_trial:.4f})")
print("This implies that the required geometric ratio (h/a) for the design is the same as the trial's ratio.")
# Set h_design / a_design = h_trial / a_plate
h_over_a_design = h_over_a_trial
h_design_m = h_over_a_design * a_design
h_design_mm = h_design_m * 1000
print(f"\nFinal Calculation:")
print(f"Required pavement thickness (h_design) = (h/a ratio) * a_design")
print(f"h_design = {h_over_a_design:.3f} * {a_design:.4f} m = {h_design_m:.4f} m")
print(f"h_design = {h_over_a_design:.3f} * {a_design*1000:.2f} mm = {h_design_mm:.2f} mm")

print(f"\nThe required pavement thickness is {h_design_mm:.2f} mm.")

# Final answer in the requested format
print(f"\n<<<{h_design_mm:.2f}>>>")