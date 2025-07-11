# Given parameters
eta_transformer = 0.98  # Transformer efficiency
loss_line_percent = 2.0   # Transmission line resistive losses in percent

# Calculate the efficiency of the transmission line
eta_line = 1 - (loss_line_percent / 100)

# Calculate the overall coefficient for delivered power
# This represents the combined effect of transformer and line efficiencies.
power_delivery_coeff = eta_transformer * eta_line

# Line impedance values
Z_line_R = 0.08  # Resistance in Ohms
Z_line_X = 0.16  # Reactance in Ohms

# The coefficient for reactive power compensation as given in the correct option
q_comp_coeff = 0.979

# Print the final equations for the user
print("Based on the analysis, the correct equations are:")
print("-" * 50)

# Print the equation for Total Real Power Delivered
# The format string uses the calculated coefficient.
print(f"P_delivered = {power_delivery_coeff:.4f} * (P_wind + P_pv)")

# Print the equation for Voltage Drop
# The format string uses the derived and given coefficients and impedance values.
print(f"V_drop = (({power_delivery_coeff:.4f} * (P_wind + P_pv}) + j * Q_comp * {q_comp_coeff}) / V_nominal) * ({Z_line_R} + j{Z_line_X})")
print("-" * 50)
