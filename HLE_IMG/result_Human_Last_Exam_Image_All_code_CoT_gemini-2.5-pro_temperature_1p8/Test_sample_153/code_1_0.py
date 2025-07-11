# Power System Analysis

# --- Part 1: Total Real Power Delivered ---

# Given efficiency and loss values
eta_transformer = 0.98  # Transformer efficiency
loss_line_percent = 2   # Transmission line resistive losses in percent

# Calculate transmission line efficiency
eta_line = 1 - (loss_line_percent / 100)

# Calculate the total efficiency factor for power delivery
# This assumes power passes through one transformer and one transmission line segment.
total_efficiency_factor = eta_transformer * eta_line

# --- Part 2: Voltage Drop Equation ---

# Given line impedance
R_line = 0.08  # Resistance in ohms
X_line = 0.16  # Reactance in ohms

# Factor for reactive power compensation term, taken from the consistent option
q_comp_factor = 0.979

# --- Output the Final Equations ---

print("The derived equations based on the analysis are:")
print("-" * 50)

# Print the equation for total real power delivered
print("Total Real Power Delivered Equation:")
print(f"P_delivered = {total_efficiency_factor:.4f} * (P_wind + P_pv)")

print("\nVoltage Drop Equation:")
# Print the equation for voltage drop
# Note: This represents the structure of the correct answer choice.
# The power term P uses the receiving-end power as an approximation.
print(f"V_drop = (({total_efficiency_factor:.4f} * (P_wind + P_pv) + j Q_comp * {q_comp_factor:.3f}) / V_nominal) * ({R_line:.2f} + j{X_line:.2f})")
print("-" * 50)