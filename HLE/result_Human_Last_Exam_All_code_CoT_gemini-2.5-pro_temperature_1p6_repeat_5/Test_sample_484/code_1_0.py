# This script prints the derived mathematical expression for the change in synaptic efficacy.
# It includes the definition of the new variables introduced in the simplification.

# Define symbolic representations for the variables in the equation.
tau_w = "\\tau_w"
w_i_dot = "\\frac{dw_i}{dt}"
u_i = "u_i"
v_i = "v_i"
beta = "\\beta"
rho = "\\rho"

# Print the final derived equation.
# The equation shows how the synaptic weight w_i changes over time.
# The change depends on a postsynaptic factor (u_i) and a presynaptic factor (v_i).
print("The derived expression for the dynamics of the synaptic efficacy w_i is:")
print(f"{tau_w} {w_i_dot} = {u_i} * ({beta} - {rho} / (1 + {v_i}))")

# Provide the definitions for the variables used in the final equation.
print("\nWhere the variables and constants are defined as:")
print(f"w_i: The synaptic efficacy (strength) of the i-th synapse.")
print(f"{v_i}: The presynaptic accumulator, representing the steady-state level of presynaptic MMP9, M_i.")
print(f"{u_i}: The postsynaptic accumulator, representing the steady-state level of shared postsynaptic calcium, Y.")
print(f"{beta}: The parameter representing the strength of LTP induced by BDNF.")
print(f"{rho}: A positive constant defined as rho = ({beta} - \\alpha) * (1 - \\eta), which combines the relative strengths of LTP and LTD with the baseline proBDNF level.")