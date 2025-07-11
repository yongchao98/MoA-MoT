# Define the parameters based on the problem description and analysis.

# Transformer efficiency
eta_transformer = 0.98

# Transmission line efficiency (100% - 2% resistive loss)
eta_line = 1 - 0.02

# Line impedance components
R_line = 0.08
X_line = 0.16

# Factor associated with reactive power compensation in the given formula
q_factor = 0.979

# --- Step 1: Calculate the total efficiency for real power delivery ---
# The total efficiency is the product of the transformer and line efficiencies.
total_efficiency = eta_transformer * eta_line

# --- Step 2: Construct and print the final equations from the chosen answer ---

# Equation for Total Real Power Delivered (P_delivered)
# This equation represents the total power from the sources multiplied by the overall efficiency.
p_delivered_eq = f"P_delivered = {total_efficiency:.4f} * (P_wind + P_pv)"

# Equation for Voltage Drop (V_drop)
# This equation calculates the voltage drop based on the delivered power, reactive compensation,
# and line impedance.
v_drop_eq = (
    f"V_drop = (({total_efficiency:.4f} * (P_wind + P_pv) + j * Q_comp * {q_factor}) "
    f"/ V_nominal) * ({R_line} + j{X_line})"
)

# Print the final equations
print("The derived equations are:")
print(p_delivered_eq)
print(v_drop_eq)
