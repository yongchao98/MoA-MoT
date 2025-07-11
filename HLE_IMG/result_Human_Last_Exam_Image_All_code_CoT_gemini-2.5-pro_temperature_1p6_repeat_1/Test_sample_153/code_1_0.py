# These are the variables representing the coefficients and values from the problem analysis.

# Efficiency coefficient for total real power delivered after three stages of loss (0.98^3 is approx 0.9404)
p_delivered_coeff = 0.9404

# Efficiency coefficient for power entering the final line section, after two stages of loss (0.98^2)
p_input_coeff_for_vdrop = 0.9604

# Coefficient associated with reactive power compensation in the given formula
q_comp_coeff = 0.979

# Given line impedance values
line_R = 0.08
line_X = 0.16

# Print the final equations corresponding to the most logical answer choice (E).
# Each number from the final equations is represented by a variable above.

print("P_delivered = {} * (P_wind + P_pv)".format(p_delivered_coeff))
print("V_drop = (({} * (P_wind + P_pv) + j * Q_comp * {}) / V_nominal) * ({} + j*{})".format(p_input_coeff_for_vdrop, q_comp_coeff, line_R, line_X))