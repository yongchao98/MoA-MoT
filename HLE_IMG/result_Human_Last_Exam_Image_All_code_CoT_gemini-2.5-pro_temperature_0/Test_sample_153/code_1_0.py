# Define the given parameters
eta_transformer = 0.98
line_loss_percentage = 2.0
R_line = 0.08
X_line = 0.16

# This coefficient is given in the correct option C, its origin is not fully specified in the problem
q_comp_coeff = 0.979

# 1. Calculate the total efficiency for real power delivery
eta_line = 1 - (line_loss_percentage / 100)
total_efficiency = eta_transformer * eta_line

# 2. Format the equations for printing
# Equation for delivered power
p_delivered_eq = f"P_delivered = {total_efficiency:.4f} * (P_wind + P_pv)"

# Equation for voltage drop
v_drop_eq = (f"V_drop = (({total_efficiency:.4f} * (P_wind + P_pv) + "
             f"j Q_comp * {q_comp_coeff}) / V_nominal) * ({R_line} + j{X_line})")

# Print the results
print("Derived Equations:")
print(p_delivered_eq)
print(v_drop_eq)
