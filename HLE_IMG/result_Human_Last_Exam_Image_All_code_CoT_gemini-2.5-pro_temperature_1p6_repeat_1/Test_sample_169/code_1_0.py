import cmath
import math

def solve_wind_compensation():
    """
    Solves for the optimal reactive power injection to restore voltage in a
    wind energy system after a single line-to-ground fault.
    """
    # 1. Define given parameters
    S_base_MVA = 100.0  # MVA
    V_base_kV = 0.575  # kV (Nominal voltage at Bus-W)

    Z_WF_ohm = 0.01 + 0.05j  # Complex impedance in Ohms
    R_F_ohm = 0.1  # Fault resistance in Ohms

    harmonic_loss_increase = 0.06  # 6% increase
    min_power_factor = 0.95

    # 2. Calculate base values for the per-unit system
    S_base_VA = S_base_MVA * 1e6
    V_base_V = V_base_kV * 1e3
    Z_base_ohm = V_base_V**2 / S_base_VA

    # 3. Convert system parameters to per-unit
    Z_WF_pu = Z_WF_ohm / Z_base_ohm
    R_F_pu = R_F_ohm / Z_base_ohm

    # 4. Calculate power required to maintain V_W = 1.0 p.u.
    # We assume Bus-W voltage is the reference V_W = 1.0 + j0.0 p.u.
    V_W_pu = 1.0 + 0.0j

    # The system is a voltage divider: V_F = V_W * R_F / (Z_WF + R_F)
    Z_total_pu = Z_WF_pu + R_F_pu
    
    # We can directly calculate the current into the fault, assuming V_G is disconnected
    # I_WF = V_W / (Z_WF + R_F) -> this is the total current from V_W to ground
    I_WF_pu = V_W_pu / Z_total_pu
    
    # Calculate apparent power drawn from Bus-W: S = V * I*
    S_W_pu = V_W_pu * I_WF_pu.conjugate()

    P_W_pu = S_W_pu.real
    Q_load_pu = S_W_pu.imag

    # 5. Account for harmonic losses
    # The real power demand increases by the harmonic loss factor
    P_gen_pu = P_W_pu * (1 + harmonic_loss_increase)

    # 6. Solve the optimization problem
    # Objective: Minimize Q_comp
    # Constraint: Q_wind + Q_comp = Q_load_pu
    # To minimize Q_comp, we must maximize Q_wind
    # The maximum Q_wind is limited by the power factor of the generator
    
    # Calculate max angle from power factor: arccos(0.95)
    max_angle_rad = math.acos(min_power_factor)
    
    # Calculate max Q_wind based on P_gen and PF limit: Q = P * tan(angle)
    Q_wind_max_pu = P_gen_pu * math.tan(max_angle_rad)

    # The optimal Q_comp is the deficit between what the load needs and what the wind farm can provide
    Q_opt_pu = Q_load_pu - Q_wind_max_pu

    # 7. Convert the optimal p.u. value back to MVAR
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    
    # Print the results and the final equation
    print("--- System Parameters (Per-Unit) ---")
    print(f"Base Impedance (Z_base): {Z_base_ohm:.6f} Ohms")
    print(f"Line Impedance (Z_WF): {Z_WF_pu.real:.4f} + j{Z_WF_pu.imag:.4f} p.u.")
    print(f"Fault Resistance (R_F): {R_F_pu:.4f} p.u.")
    print("")
    print("--- Power Requirement Calculation ---")
    print(f"Fundamental Real Power Load (P_load): {P_W_pu:.5f} p.u.")
    print(f"Total Real Power to Generate (P_gen with harmonic loss): {P_gen_pu:.5f} p.u.")
    print(f"Reactive Power Load (Q_load): {Q_load_pu:.5f} p.u.")
    print("")
    print("--- Optimization Solution ---")
    print(f"Max reactive power from wind farm at 0.95 PF (Q_wind_max): {Q_wind_max_pu:.5f} p.u.")
    
    print("\nThe optimal reactive power injection (Q_opt) is found by:")
    print("Q_opt = Q_load - Q_wind_max")
    # Print the final equation with numbers
    print(f"Q_opt = {Q_load_pu:.5f} p.u. - {Q_wind_max_pu:.5f} p.u. = {Q_opt_pu:.5f} p.u.")

    print("\nConverting the result to MVAR:")
    print(f"Q_opt_MVAR = Q_opt_p.u. * S_base_MVA")
    print(f"Q_opt_MVAR = {Q_opt_pu:.5f} * {S_base_MVA} MVA = {Q_opt_MVAR:.4f} MVAR")

    print("\n--- Final Answer ---")
    print(f"The optimal reactive power injection required from the compensating device is {Q_opt_MVAR:.4f} MVAR.")
    
    # Return the final numeric answer for the platform
    return Q_opt_MVAR

# Execute the function and capture the result
final_answer = solve_wind_compensation()
# Present the final answer in the required format
# The value is formatted to one decimal place as requested in the example format.
print(f"\n<<<{final_answer:.1f}>>>")
