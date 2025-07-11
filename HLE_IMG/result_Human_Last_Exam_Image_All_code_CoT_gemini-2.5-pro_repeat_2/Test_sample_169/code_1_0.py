import numpy as np
from scipy.optimize import fsolve

# --- Step 1: Define System Parameters in Per Unit (p.u.) ---
# Base values
S_base = 100e6  # VA
# The problem statement is ambiguous about the voltage base.
# The most plausible interpretation for a textbook problem is that
# the impedance and resistance values are already given in p.u.
# Z_WF = 0.01 + j0.05 p.u.
# R_F = 0.1 p.u.
R_WF = 0.01  # p.u.
X_WF = 0.05  # p.u.
R_F = 0.1    # p.u.

# System conditions and targets
# The target voltage "0.575 kV" is interpreted as 0.575 p.u. for a feasible problem setup.
V_W_target = 0.575  # p.u. Target voltage at Bus-W
# Assume a constant active power injection from the wind farm. A typical value for
# a large farm would be near rated power. We assume 80 MW.
P_W = 80 / (S_base / 1e6) # p.u. Active power from wind farm
Q_max_MVAR = 10
Q_max = Q_max_MVAR / (S_base / 1e6) # p.u. Maximum reactive power injection
pf_min = 0.95  # Minimum required power factor (lagging)
harmonic_loss_factor = 1.06 # 6% increase in total system losses

# --- Step 2: Formulate the System Model ---
# We use a simplified model representing the grid as a Thévenin equivalent
# with voltage V_G = 1.0 p.u. The impedance from Bus-W to this source
# includes the line and the fault.
# Z_eq = Z_WF + R_F = (R_WF + R_F) + jX_WF
R_base = R_WF + R_F
# The harmonic losses increase the resistive component of the impedance.
R_eff = R_base * harmonic_loss_factor
X_eff = X_WF

# The power flow equations relate the voltages and powers:
# V_G = V_W - I_W * Z_eq  (assuming current flows from W to G)
# This can be expanded into two real equations using V_W = V * e^(j*delta)
# To solve for Q_comp that achieves V = V_W_target, we eliminate delta
# using cos^2 + sin^2 = 1, leading to the equation below.

def voltage_equation(Q_comp, V, P, R, X):
    """
    Represents the power flow equation that needs to be solved for Q_comp.
    Returns the difference that should be zero for a valid solution.
    (V^2 - P*R - Q*X)^2 + (P*X - Q*R)^2 - V^2 = 0
    """
    term1 = (V**2 - P * R - Q_comp * X)**2
    term2 = (P * X - Q_comp * R)**2
    return term1 + term2 - V**2

# --- Step 3: Solve for the Required Reactive Power ---
# Use a numerical solver to find the Q_comp that satisfies the equation for V = V_W_target.
initial_guess = 0.1
Q_required_pu = fsolve(voltage_equation, initial_guess, args=(V_W_target, P_W, R_eff, X_eff))[0]
Q_required_MVAR = Q_required_pu * (S_base / 1e6)

# --- Step 4: Analyze the Solution and Constraints ---
# Check if the required Q is within the device's capability
is_feasible_q = 0 <= Q_required_pu <= Q_max

# Check the power factor constraint
# pf = P_W / sqrt(P_W^2 + Q_comp^2) >= pf_min -> Q_comp^2 <= P_W^2 * (1/pf_min^2 - 1)
Q_pf_max_pu = P_W * np.sqrt(1/pf_min**2 - 1)
is_feasible_pf = Q_required_pu <= Q_pf_max_pu

# --- Step 5: Determine the Optimal Q and Print Results ---
print("--- System Parameters ---")
print(f"Target Voltage |V_W|: {V_W_target:.3f} p.u.")
print(f"Wind Farm Active Power P_W: {P_W:.2f} p.u. ({P_W * 100} MW)")
print(f"Effective Line Resistance (with harmonic losses) R_eff: {R_eff:.4f} p.u.")
print(f"Effective Line Reactance X_eff: {X_eff:.4f} p.u.")
print(f"Max Reactive Power Q_max: {Q_max:.2f} p.u. ({Q_max_MVAR} MVAR)")
print("-" * 25)

print("--- Problem Formulation ---")
print("The relationship between voltage |V_W|, active power P_W, and reactive power Q_comp is:")
print("(|V_W|^2 - P_W*R_eff - Q_comp*X_eff)^2 + (P_W*X_eff - Q_comp*R_eff)^2 - |V_W|^2 = 0")
print("\nSubstituting the values, we need to solve for Q_comp in the equation:")
term1_const = V_W_target**2 - P_W*R_eff
term1_q_coeff = -X_eff
term2_const = P_W*X_eff
term2_q_coeff = -R_eff
rhs = V_W_target**2
print(f"({term1_const:.4f} + ({term1_q_coeff:.4f})*Q_comp)^2 + ({term2_const:.4f} + ({term2_q_coeff:.4f})*Q_comp)^2 = {rhs:.4f}")
print("-" * 25)

print("--- Analysis Results ---")
print(f"Required Q to reach {V_W_target:.3f} p.u.: {Q_required_pu:.4f} p.u. ({Q_required_MVAR:.2f} MVAR)")

if Q_required_pu > Q_max:
    print(f"This required Q ({Q_required_pu:.4f} p.u.) exceeds the compensator's maximum capacity of {Q_max:.2f} p.u.")
    print("Therefore, it is not possible to restore the voltage to the target of 0.575 p.u.")
    print("The optimization problem is to minimize Q_comp while trying to meet the voltage target.")
    print("Since the target is unreachable, the best strategy to maximize voltage support is to inject the maximum available reactive power.")
    Q_opt_pu = Q_max
    print(f"\nCONCLUSION: The optimal reactive power is the maximum available from the device.")
else:
    # This case is unlikely based on the problem numbers but included for completeness.
    if Q_required_pu <= Q_pf_max_pu:
        Q_opt_pu = Q_required_pu
        print("This required Q is feasible and satisfies all constraints.")
        print(f"\nCONCLUSION: The optimal reactive power is the minimum required.")
    else:
        Q_opt_pu = Q_pf_max_pu
        print("This required Q violates the power factor constraint.")
        print(f"Maximum Q allowed by PF constraint: {Q_pf_max_pu:.4f} p.u.")
        print(f"\nCONCLUSION: The optimal reactive power is limited by the power factor requirement.")


Q_opt_MVAR = Q_opt_pu * (S_base / 1e6)
print("-" * 25)
print(f"Final Optimal Reactive Power Injection Q_opt = {Q_opt_pu:.4f} p.u. = {Q_opt_MVAR:.2f} MVAR")
<<<10.0>>>