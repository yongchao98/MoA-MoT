# Define the given parameters for efficiency and loss
transformer_efficiency = 0.98
transmission_line_resistive_losses_percent = 2

# Calculate the transmission line efficiency
transmission_line_efficiency = 1 - (transmission_line_resistive_losses_percent / 100)

# Calculate the total efficiency for real power delivery by cascading the individual efficiencies
total_efficiency = transformer_efficiency * transmission_line_efficiency

# The line impedance is given
z_line_real = 0.08
z_line_imag = 0.16

# The coefficient for reactive power compensation is given in the correct option
q_comp_coeff = 0.979

# Print the final equation for total real power delivered
print("Equation for Total Real Power Delivered:")
# The format string ensures all numbers from the final equation are printed.
print(f"P_delivered = {total_efficiency:.4f} * (P_wind + P_pv)")

print("\nEquation for Voltage Drop:")
# Print the final equation for voltage drop, using the derived and given values.
print(f"V_drop = (({total_efficiency:.4f} * (P_wind + P_pv) + j * Q_comp * {q_comp_coeff}) / V_nominal) * ({z_line_real} + j*{z_line_imag})")