# Define the variables symbolically for printing the equations.
P_wind = "P_wind"
P_pv = "P_pv"
Q_comp = "Q_comp"
V_nominal = "V_nominal"

# 1. Calculate the efficiency coefficient for delivered power
transformer_efficiency = 0.98
line_loss_efficiency = 1 - 0.02
total_efficiency = transformer_efficiency * line_loss_efficiency

# 2. Define the line impedance and the reactive power coefficient from the correct option
line_impedance_R = 0.08
line_impedance_X = 0.16
q_comp_coeff = 0.979

# 3. Print the first equation for delivered power
# The format ensures all numbers from the final equation are printed.
print(f"P_delivered = {total_efficiency:.4f} * ({P_wind} + {P_pv})")

# 4. Print the second equation for voltage drop
# The format ensures all numbers from the final equation are printed.
print(f"V_drop = (({total_efficiency:.4f} * ({P_wind} + {P_pv}) + j * {Q_comp} * {q_comp_coeff}) / {V_nominal}) * ({line_impedance_R} + j{line_impedance_X})")
