import cmath

# Define the constants based on the problem description and analysis.
transformer_efficiency = 0.98
line_loss_percentage = 0.02
line_efficiency = 1 - line_loss_percentage

# Calculate the total efficiency for real power delivery.
total_efficiency = transformer_efficiency * line_efficiency

# Line impedance values from the problem statement.
line_resistance = 0.08
line_reactance = 0.16

# Reactive power coefficient from the selected option.
q_comp_factor = 0.979

# Generate and print the equations from the correct answer choice.
print("# Equation for Total Real Power Delivered:")
# The format f"{value:.4f}" ensures the number is printed with 4 decimal places.
equation1 = f"    P_delivered = {total_efficiency:.4f} * (P_wind + P_pv)"
print(equation1)
print("\n# Equation for Voltage Drop:")
equation2 = f"    V_drop = (({total_efficiency:.4f} * (P_wind + P_pv) + j Q_comp * {q_comp_factor}))/V_nominal * ({line_resistance} + j{line_reactance})"
print(equation2)
