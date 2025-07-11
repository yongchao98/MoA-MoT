import math

# --- Given Parameters ---
V_W_target = 0.575  # Target voltage at Bus-W (p.u.)
Z_WF_R = 0.01      # Transmission line resistance (p.u.)
Z_WF_X = 0.05      # Transmission line reactance (p.u.)
R_F = 0.1          # Fault resistance (p.u.)
S_base = 100.0     # Base power (MVA)
Q_max = 10.0       # Maximum reactive power of compensator (MVAR)
PF_min = 0.95      # Minimum power factor (lagging)
loss_factor = 1.06 # Factor for harmonic losses (1 + 6%)

def solve_optimal_q():
    """
    Calculates the optimal reactive power injection to restore bus voltage.
    """
    print("Step 1: Define system parameters in per-unit.")
    Q_max_pu = Q_max / S_base
    print(f"  - Target Voltage |V_W|: {V_W_target} p.u.")
    print(f"  - Line Impedance Z_WF: {Z_WF_R} + j{Z_WF_X} p.u.")
    print(f"  - Fault Resistance R_F: {R_F} p.u.")
    print(f"  - Compensator Q_max: {Q_max_pu:.2f} p.u. ({Q_max} MVAR)")
    print("-" * 30)

    print("Step 2: Calculate the power drawn by the faulted line to maintain target voltage.")
    # Total impedance from Bus-W to ground through the fault
    R_total = Z_WF_R + R_F
    # Denominator for power calculation: |Z_total|^2 = R_total^2 + X_WF^2
    Z_total_sq = R_total**2 + Z_WF_X**2
    
    # Power drawn by the faulted line: S_line = V^2 / Z*
    # P_line = Re{V^2 / (R_total - jX_WF)} = V^2 * R_total / |Z_total|^2
    P_line = (V_W_target**2 * R_total) / Z_total_sq
    # Q_line = Im{V^2 / (R_total - jX_WF)} = V^2 * X_WF / |Z_total|^2
    Q_line = (V_W_target**2 * Z_WF_X) / Z_total_sq
    print(f"  - Active power required by line (P_line): {P_line:.4f} p.u.")
    print(f"  - Reactive power required by line (Q_line): {Q_line:.4f} p.u.")
    print("-" * 30)
    
    print("Step 3: Account for harmonic losses.")
    # Total active power generation must cover line power + harmonic losses
    P_gen = P_line * loss_factor
    print(f"  - Total active power generation needed (P_gen): {P_line:.4f} * {loss_factor} = {P_gen:.4f} p.u.")
    print("-" * 30)
    
    print("Step 4: Optimize Q_comp by maximizing Q_gen from the wind farm.")
    # To minimize Q_comp, we must maximize Q_gen subject to the PF constraint.
    # PF = P_gen / |S_gen| >= PF_min  => |S_gen| <= P_gen / PF_min
    # Q_gen = sqrt(|S_gen|^2 - P_gen^2)
    # The max Q_gen occurs at the minimum PF.
    # tan(phi_max) = sqrt(1/PF_min^2 - 1)
    tan_phi_max = math.sqrt(1 / PF_min**2 - 1)
    Q_gen_max = P_gen * tan_phi_max
    print(f"  - Generator PF constrained to >= {PF_min} lagging.")
    print(f"  - Maximum reactive power from generator (Q_gen): {Q_gen_max:.4f} p.u.")
    print("-" * 30)

    print("Step 5: Solve for the optimal reactive power injection (Q_opt).")
    # Reactive power balance at Bus-W: Q_line = Q_gen + Q_comp
    Q_opt_pu = Q_line - Q_gen_max
    Q_opt_MVAR = Q_opt_pu * S_base
    
    print("The final reactive power balance equation is:")
    print(f"  Q_line = Q_gen + Q_comp")
    print(f"  {Q_line:.4f} p.u. = {Q_gen_max:.4f} p.u. + Q_opt")
    print(f"\nThis gives the optimal reactive power injection:")
    print(f"  Q_opt = {Q_line:.4f} - {Q_gen_max:.4f} = {Q_opt_pu:.4f} p.u.")
    print(f"  Q_opt = {Q_opt_MVAR:.2f} MVAR")
    print("-" * 30)

    print("Step 6: Final feasibility check.")
    if Q_opt_MVAR > Q_max:
        print(f"  WARNING: The required optimal injection ({Q_opt_MVAR:.2f} MVAR) exceeds the device's maximum capacity of {Q_max} MVAR.")
        print("  The target voltage of 0.575 p.u. is not achievable with the given constraints.")
    else:
        print(f"  The required injection ({Q_opt_MVAR:.2f} MVAR) is within the device's maximum capacity of {Q_max} MVAR.")

    return Q_opt_MVAR

# --- Execute the solution ---
if __name__ == "__main__":
    optimal_q = solve_optimal_q()
    # The final numerical answer is extracted for the platform.
    # Note: The calculation shows this is the *required* Q to meet the goal, 
    # even if it exceeds the device's physical limit.
    final_answer = optimal_q
    # print(f"\nFinal Answer: {final_answer:.2f}") # This is for display, not the final return
    # print(f"<<<{final_answer:.2f}>>>")


solve_optimal_q()