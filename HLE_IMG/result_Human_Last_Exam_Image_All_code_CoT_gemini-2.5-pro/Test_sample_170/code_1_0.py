import cmath
import math

# Step 1: Define system parameters in per-unit (pu)
S_base_MVA = 100.0  # MVA
V_base_kV = 220.0  # kV
V_th_pu = 1.0  # Thevenin voltage (grid source) in pu
R_th_pu = 0.02  # Thevenin resistance (system impedance) in pu
X_th_pu = 0.10  # Thevenin reactance (system impedance) in pu
Z_th_pu = complex(R_th_pu, X_th_pu)
V_b_sag_pu = 0.85  # Voltage at Bus B during fault in pu
V_b_target_pu = 1.0  # Target voltage at Bus B in pu
harmonic_loss_increase = 0.04  # 4% increase in losses

# Step 2: Analyze the initial state to find the equivalent load PL (assuming QL=0)
# The voltage-power relationship is a quadratic equation of the form:
# a*PL^2 + b*PL + c = 0
# where:
# a = R_th^2 + X_th^2
# b = 2 * R_th * V_b^2
# c = V_b^4 - V_b^2 * V_th^2
print("Step 1: Calculate the equivalent real power load (PL) causing the voltage sag to 0.85 pu.")
a_pl = R_th_pu**2 + X_th_pu**2
b_pl = 2 * R_th_pu * V_b_sag_pu**2
c_pl = V_b_sag_pu**4 - (V_b_sag_pu**2 * V_th_pu**2)

print(f"Solving the quadratic equation for PL: {a_pl:.4f}*PL^2 + {b_pl:.4f}*PL + {c_pl:.4f} = 0")

# Solve quadratic equation for PL
discriminant_pl = b_pl**2 - 4 * a_pl * c_pl
pl_sol1 = (-b_pl + math.sqrt(discriminant_pl)) / (2 * a_pl)
pl_sol2 = (-b_pl - math.sqrt(discriminant_pl)) / (2 * a_pl)

# Choose the physically meaningful positive power solution
P_L_pu = pl_sol1 if pl_sol1 > 0 else pl_sol2
print(f"The equivalent real power load is PL = {P_L_pu:.4f} pu.\n")


# Step 3: Calculate the required reactive power injection (Q_opt)
# In the final state, V_b = 1.0 pu. The net reactive power at the bus is Q_net = QL - Q_opt.
# With our assumption QL=0, Q_net = -Q_opt.
# This gives a quadratic equation for Q_net (and thus for Q_opt):
# a*Q_net^2 + b*Q_net + c = 0
# where:
# a = R_th^2 + X_th^2
# b = 2 * X_th * V_b^2
# c = V_b^4 - V_b^2*V_th^2 + 2*P_L*R_th*V_b^2 + P_L^2*(R_th^2+X_th^2)
# Simplified form:
# a_q*Q_net^2 + b_q*Q_net + c_q = 0
print("Step 2: Calculate the optimal reactive power (Q_opt) to restore voltage to 1.0 pu.")
a_q = R_th_pu**2 + X_th_pu**2
b_q = 2 * X_th_pu * V_b_target_pu**2
c_q = (2 * P_L_pu * R_th_pu * V_b_target_pu**2) + (P_L_pu**2 * (R_th_pu**2 + X_th_pu**2))

print(f"Solving the quadratic equation for Q_net: {a_q:.4f}*Q_net^2 + {b_q:.4f}*Q_net + {c_q:.4f} = 0")

# Solve quadratic equation for Q_net
discriminant_q = b_q**2 - 4 * a_q * c_q
q_net_sol1 = (-b_q + math.sqrt(discriminant_q)) / (2 * a_q)
q_net_sol2 = (-b_q - math.sqrt(discriminant_q)) / (2 * a_q)

# Q_opt = -Q_net. We want the minimum positive Q_opt.
Q_opt_sol1_pu = -q_net_sol1
Q_opt_sol2_pu = -q_net_sol2

Q_opt_pu = min(Q_opt_sol1_pu, Q_opt_sol2_pu)
Q_opt_MVAR = Q_opt_pu * S_base_MVA
print(f"The minimum required reactive power injection is Q_opt = {Q_opt_pu:.4f} pu, or {Q_opt_MVAR:.2f} MVAR.\n")

# Step 4 & 5: Calculate system losses with harmonic effects
print("Step 3: Calculate the total system power losses.")
# Total complex power drawn from the bus in the final state
S_net_pu = complex(P_L_pu, -Q_opt_pu)
# Current from the source
# S = V*I_conj => I_conj = S/V => I = (S/V)_conj = S_conj/V_conj
# Since V_B is the reference, lets use I = (V_th - V_b)/Z_th
# To calculate angle delta of V_b:
# S_net = V_b * conj((V_th - V_b)/Z_th)
# We can find current magnitude from power: |I| = |S_net| / |V_b_target|
I_mag_pu = abs(S_net_pu) / V_b_target_pu
P_loss_base_pu = I_mag_pu**2 * R_th_pu
P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
P_loss_total_MW = P_loss_total_pu * S_base_MVA

print(f"The current magnitude is |I| = sqrt({P_L_pu:.4f}^2 + {Q_opt_pu:.4f}^2) / {V_b_target_pu:.1f} = {I_mag_pu:.4f} pu.")
print(f"Base system losses = I^2 * R = {I_mag_pu:.4f}^2 * {R_th_pu:.2f} = {P_loss_base_pu:.4f} pu.")
print(f"Total system losses with 4% increase = {P_loss_base_pu:.4f} * (1 + {harmonic_loss_increase:.2f}) = {P_loss_total_pu:.4f} pu, or {P_loss_total_MW:.2f} MW.\n")

# Step 6: Final summary
print("--- Summary of Results ---")
print(f"Optimal reactive power injection from MMCC STATCOM (Q_opt): {Q_opt_MVAR:.2f} MVAR")
print(f"Total system real power losses: {P_loss_total_MW:.2f} MW")
print("\nNote: The calculated reactive power injection of 126.45 MVAR exceeds the STATCOM's maximum capacity of 50 MVAR. "
      "Also, the calculation assumes an ideal STATCOM (zero real power consumption), which would violate the power factor constraint. "
      "This indicates that restoring the voltage to 1.0 pu is not feasible with the given STATCOM constraints under this fault condition.")
<<<126.45>>>