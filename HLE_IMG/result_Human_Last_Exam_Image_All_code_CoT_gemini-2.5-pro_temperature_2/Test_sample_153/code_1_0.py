# Define the coefficients based on the problem description and analysis
transformer_efficiency = 0.98
line_loss_percentage = 0.02
line_efficiency = 1 - line_loss_percentage

# The problem states the line impedance
r_line = 0.08
x_line = 0.16

# This coefficient is given in the correct option's formula
q_comp_coeff = 0.979

# Calculate the overall power delivery coefficient by multiplying the efficiencies
# P_delivered_coeff = transformer_efficiency * line_efficiency
# This corresponds to the calculation 0.98 * 0.98
p_delivered_coeff = 0.9604

# --- Output the Final Equations ---

# Print the equation for the Total Real Power Delivered
print("Total Real Power Delivered Equation:")
print(f"P_delivered = {p_delivered_coeff} * (P_wind + P_pv)")

print("\nVoltage Drop Equation:")
# Print the equation for the Voltage Drop
# Note: In a real scenario, the complex voltage drop is often calculated as I * Z.
# The formula provided in the options is a specific form we must adhere to.
# V_drop = (S / V_nominal) * Z, where S is the complex power.
# The real power P in S is taken as the P_delivered calculated above.
# The reactive power Q is given as a function of Q_comp.
print(f"V_drop = ({p_delivered_coeff} * (P_wind + P_pv) + j Q_comp * {q_comp_coeff}) / V_nominal * ({r_line} + j{x_line})")