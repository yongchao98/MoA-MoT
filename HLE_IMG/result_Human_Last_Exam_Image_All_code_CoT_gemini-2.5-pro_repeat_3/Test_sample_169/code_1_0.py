import math

def solve_wind_system_optimization():
    """
    Formulates and solves the nonlinear optimization problem to find the
    optimal reactive power injection from a fault compensating device.
    """
    # --- Given Parameters ---
    V_W_target_pu = 0.575      # Target voltage at Bus-W in p.u.
    Z_WF_pu = 0.01 + 0.05j     # Transmission impedance in p.u.
    R_F_pu = 0.1               # Fault resistance in p.u.
    Q_max_MVAR = 10            # Max reactive power of compensator in MVAR
    S_base_MVA = 100           # Base power in MVA
    PF_min_lagging = 0.95      # Minimum power factor for the generator
    harmonic_loss_factor = 0.06 # 6% of total losses are from harmonics

    # --- Step 1: Model the faulted load in per-unit ---
    # The load is the series impedance of the line segment and the fault resistance.
    R_pu = Z_WF_pu.real + R_F_pu
    X_pu = Z_WF_pu.imag
    Z_pu_mag = math.sqrt(R_pu**2 + X_pu**2)

    # --- Step 2: Formulate and solve for total power requirements at Bus-W ---
    # Total apparent power squared flowing from Bus-W: |S_W|^2 = (|V_W|^2 / |Z|)^2
    S_W_mag_sq = (V_W_target_pu**2 / Z_pu_mag)**2

    # Relate total active power (P_W) and reactive power (Q_W) considering losses.
    # Total active loss P_W = Fundamental Loss P_L + Harmonic Loss P_H
    # P_H = harmonic_loss_factor * P_W, so P_L = (1 - harmonic_loss_factor) * P_W
    # Since P_L = |I|^2 * R and Q_W = |I|^2 * X, we have Q_W / X = P_L / R.
    # This gives: Q_W = (P_L / R) * X = ((1 - harmonic_loss_factor) * P_W / R_pu) * X_pu
    k = (1 - harmonic_loss_factor) * (X_pu / R_pu) # Q_W = k * P_W

    # Solve the system: P_W^2 + Q_W^2 = |S_W|^2
    # P_W^2 + (k * P_W)^2 = S_W_mag_sq
    # P_W^2 * (1 + k^2) = S_W_mag_sq
    P_W_sq = S_W_mag_sq / (1 + k**2)
    P_W_pu = math.sqrt(P_W_sq)
    Q_W_pu = k * P_W_pu

    # --- Step 3: Determine the optimal reactive power compensation ---
    # The generator provides all active power: P_gen = P_W.
    P_gen_pu = P_W_pu

    # To minimize Q_comp, we maximize Q_gen subject to the power factor constraint.
    # tan(phi_max) = Q_gen_max / P_gen
    # phi_max = acos(PF_min)
    tan_phi_max = math.tan(math.acos(PF_min_lagging))
    Q_gen_max_pu = P_gen_pu * tan_phi_max

    # The optimal (minimum) compensation Q_opt is the deficit Q_W - Q_gen_max.
    Q_opt_pu = Q_W_pu - Q_gen_max_pu
    Q_opt_MVAR = Q_opt_pu * S_base_MVA

    # --- Step 4: Print the detailed solution ---
    print("--- Nonlinear Optimization Problem Solution ---")
    print(f"\nObjective: Determine the minimum reactive power injection Q_opt.")
    print("\nConstraints & Conditions:")
    print(f"1. Restore Bus-W voltage to |V_W| = {V_W_target_pu} p.u.")
    print(f"2. Generator power factor must be >= {PF_min_lagging} lagging.")
    print(f"3. Harmonic effects account for {harmonic_loss_factor*100}% of total losses.")

    print("\n--- Calculation Results ---")
    print("\n1. Total Power Required at Bus-W:")
    print(f"   - Active Power P_W   = {P_W_pu:.4f} p.u. ({P_W_pu * S_base_MVA:.2f} MW)")
    print(f"   - Reactive Power Q_W = {Q_W_pu:.4f} p.u. ({Q_W_pu * S_base_MVA:.2f} MVAR)")

    print("\n2. Generator's Maximum Reactive Power Contribution:")
    print(f"   - Active Power P_gen = {P_gen_pu:.4f} p.u.")
    print(f"   - Max Q_gen at PF={PF_min_lagging} = {Q_gen_max_pu:.4f} p.u. ({Q_gen_max_pu * S_base_MVA:.2f} MVAR)")

    print("\n3. Final Equation for Optimal Reactive Power Injection (Q_opt):")
    print("   Q_opt = Q_W - Q_gen_max")
    print(f"   In per-unit: Q_opt = {Q_W_pu:.4f} - {Q_gen_max_pu:.4f} = {Q_opt_pu:.4f} p.u.")
    print(f"   In MVAR:     Q_opt = {Q_W_pu * S_base_MVA:.4f} - {Q_gen_max_pu * S_base_MVA:.4f} = {Q_opt_MVAR:.4f} MVAR")

    print("\n--- Final Answer ---")
    print(f"The optimal reactive power injection required from the compensating device is {Q_opt_MVAR:.4f} MVAR.")
    print(f"(Note: This required injection of {Q_opt_MVAR:.2f} MVAR exceeds the device's maximum capacity of {Q_max_MVAR} MVAR, indicating the target voltage may not be achievable with the current hardware.)")

    return Q_opt_MVAR

# Execute the function to solve and print the results.
final_answer = solve_wind_system_optimization()
#<<<24.8115>>>