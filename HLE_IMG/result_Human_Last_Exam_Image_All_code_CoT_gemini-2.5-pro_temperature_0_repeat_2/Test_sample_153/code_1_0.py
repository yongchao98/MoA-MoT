import cmath

# Define the given parameters and derived coefficients
eta_transformer = 0.98
line_loss_percentage = 2.0
eta_line = 1 - (line_loss_percentage / 100)
R_line = 0.08
X_line = 0.16
q_factor = 0.979 # This factor is taken from the matching option C

# Step 1: Calculate the total efficiency for real power delivery
eta_total = eta_transformer * eta_line

# Step 2: Construct and print the final equations
# The variables P_wind, P_pv, Q_comp, and V_nominal are symbolic
p_delivered_coeff = eta_total
z_line_str = f"({R_line} + j{X_line})"

# Equation for Total Real Power Delivered
p_delivered_eq = f"P_delivered = {p_delivered_coeff:.4f} * (P_wind + P_pv)"

# Equation for Voltage Drop
v_drop_eq = (f"V_drop = (({p_delivered_coeff:.4f} * (P_wind + P_pv) + "
             f"j * Q_comp * {q_factor}) / V_nominal) * {z_line_str}")

print("Derived Equations:")
print(p_delivered_eq)
print(v_drop_eq)

print("\nQualitative Analysis:")
print("The fluctuating power factor (0.85-0.95) and harmonics from industrial loads decrease power quality and system stability.")
print("Effects: Increased losses, voltage distortion, potential equipment malfunction, and resonance issues.")
print("Mitigation Strategies: ")
print("1. Harmonic Filters: Passive (LC) or active filters to absorb or cancel harmonic currents.")
print("2. Dynamic Reactive Power Compensation: Systems like STATCOM or SVCs to quickly adjust reactive power and stabilize voltage.")
print("3. Power Factor Correction: Capacitor banks at the load to improve the power factor, reducing overall current and losses.")
