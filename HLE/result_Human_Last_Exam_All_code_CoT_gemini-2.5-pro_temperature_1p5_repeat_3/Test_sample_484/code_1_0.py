# The final expression is derived as tau_w * dw_i/dt = beta * (B_i - rho * P_i)
# We identify u_i with the BDNF level B_i and v_i with the proBDNF level P_i.

# Define the symbols for the equation
tau_w = "τ_w"
w_i_dot = "dw_i/dt"
beta = "β"
u_i = "u_i"
rho = "ρ"
v_i = "v_i"

# Print the final expression part by part
print(f"The derived expression for the change in synaptic efficacy is:")
print(f"{tau_w} * {w_i_dot} = {beta} * ({u_i} - {rho} * {v_i})")