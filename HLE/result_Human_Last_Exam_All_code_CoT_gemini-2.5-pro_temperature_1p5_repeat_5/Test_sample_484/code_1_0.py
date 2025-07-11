# This script presents the result of a steady-state analysis on the provided
# biophysical model. The derivation is outlined in the comments.

# Derivation Steps:
# 1. We start with the assumption that the variables M_i, Y, P_i, and B_i
#    reach equilibrium much faster than the synaptic weights W_i. We can therefore
#    set their time derivatives to 0 to find their steady-state values.
#    Let r_i be the average firing rate for synapse i, replacing x_i(t).
#
# 2. Solving for steady states (denoted with a '_ss' suffix):
#    - From dM_i/dt = 0: M_i_ss = phi * r_i
#    - From dY/dt = 0:   Y_ss = sum over j of (w_j * r_j)
#    - Summing the equations for dP_i/dt and dB_i/dt gives d(P_i+B_i)/dt, leading to the
#      steady-state relation: P_i_ss + B_i_ss = Y_ss.
#
# 3. We define the new variables as requested:
#    - Synaptic efficacy: w_i = W_i
#    - Presynaptic accumulator: v_i = M_i_ss
#    - Postsynaptic accumulator: u_i = Y_ss (note: u_i is the same for all synapses)
#
# 4. We substitute the steady-state relations into the weight dynamics equation:
#    tau_w * dw_i/dt = alpha * P_i + beta * B_i
#
#    Using B_i_ss = u_i - P_i_ss, we get:
#    tau_w * dw_i/dt = alpha * P_i_ss + beta * (u_i - P_i_ss)
#                   = beta * u_i + (alpha - beta) * P_i_ss
#
#    From the dP_i/dt = 0 equation, we find:
#    P_i_ss = (1 - eta) * u_i / (1 + v_i)
#
#    Substituting this back gives:
#    tau_w * dw_i/dt = beta * u_i + (alpha - beta) * (1 - eta) * u_i / (1 + v_i)
#
# 5. Simplifying the expression algebraically:
#    tau_w * dw_i/dt = u_i * [ beta * (1 + v_i) + (alpha - beta) * (1 - eta) ] / (1 + v_i)
#                   = u_i * [ beta*v_i + alpha*(1-eta) + beta*eta ] / (1 + v_i)
#                   = (beta * u_i / (1 + v_i)) * [ v_i + (alpha*(1-eta) + beta*eta) / beta ]
#
# 6. We define the constant rho to represent the plasticity threshold, which simplifies the term in the brackets to (v_i - rho):
#    rho = - (alpha*(1-eta) + beta*eta) / beta
#    This leads to the final compact expression.

print("Based on the steady-state analysis, the system can be reduced. The final expression involves a constant 'rho', defined as:")
print("\rho = - (alpha * (1 - eta) + beta * eta) / beta")

print("\nThe resulting expression for the dynamics of synaptic efficacy, tau_w * dw_i/dt, is:")
print("tau_w * dw_i/dt = (beta * u_i * (v_i - rho)) / (1 + v_i)")