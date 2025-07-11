# This script states the derived equation from the steady-state analysis.
# The goal is to express tau_w * dw_i/dt in terms of the synaptic efficacy w_i,
# the postsynaptic accumulator u_i, the presynaptic accumulator v_i, and a constant rho.

# Derived variable and constant definitions:
# v_i: The presynaptic accumulator, identified as the steady-state MMP9 level.
# u_i: The postsynaptic accumulator, identified as the shared steady-state postsynaptic calcium level.
# rho: A constant defined from the model's parameters as rho = (alpha - beta) * (1 - eta).

# The final equation combines these terms.
# Here we print each component of the final equation as a string.

tau_w = "tau_w"
dw_i_dt = "dw_i/dt"
u_i = "u_i"
beta = "beta"
rho = "rho"
one = "1"
v_i = "v_i"

print("The derived expression for the dynamics of the synaptic efficacy is:")
print(f"<<< {tau_w} * {dw_i_dt} = {u_i} * ({beta} + {rho} / ({one} + {v_i})) >>>")