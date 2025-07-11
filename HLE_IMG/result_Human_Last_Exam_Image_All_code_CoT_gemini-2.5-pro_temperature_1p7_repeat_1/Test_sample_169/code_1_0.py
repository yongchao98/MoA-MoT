import math

def solve_reactive_power_optimization():
    """
    Solves for the optimal reactive power injection from a fault compensating device.
    """
    # 1. Define given parameters in SI units
    V_target = 0.575 * 1000  # Target voltage at Bus-W in Volts
    R_WF = 0.01  # Transmission line resistance in Ohms
    X_WF = 0.05  # Transmission line reactance in Ohms
    R_F = 0.1    # Fault resistance in Ohms
    PF_min = 0.95  # Minimum lagging power factor for the generator
    Q_max = 10 * 1e6 # Maximum reactive power capacity of the device in VAR

    # 2. Calculate total impedance and required power injection based on the model
    # Model: V_W supplies a load Z_total = (R_WF + R_F) + j*X_WF
    # S_total = V_W^2 / Z_total_conjugate
    R_total = R_WF + R_F
    Z_total_sq = R_total**2 + X_WF**2
    V_target_sq = V_target**2

    # Calculate required total active power (P_g) and reactive power (Q_total)
    # P_g = Re(S_total) = Re(V^2 * Z_total / |Z_total|^2) = V^2 * R_total / |Z_total|^2
    P_g = (V_target_sq * R_total) / Z_total_sq
    # Q_total = Im(S_total) = Im(V^2 * Z_total / |Z_total|^2) = V^2 * X_WF / |Z_total|^2
    Q_total = (V_target_sq * X_WF) / Z_total_sq

    # 3. Apply the power factor constraint to find max generator reactive power (Q_g_max)
    # For a lagging PF, Q_g >= 0. PF = cos(phi), so |tan(phi)| = |Q_g/P_g|
    # To meet PF >= 0.95, |Q_g/P_g| <= tan(acos(0.95))
    tan_phi_max = math.tan(math.acos(PF_min))
    Q_g_max = P_g * tan_phi_max

    # 4. Optimize: To minimize Q_comp, maximize Q_g
    # Q_total = Q_g + Q_comp  => Q_comp = Q_total - Q_g
    # Set Q_g to its maximum possible value.
    Q_g_opt = Q_g_max
    
    # 5. Calculate the optimal compensator reactive power (Q_opt)
    Q_opt = Q_total - Q_g_opt
    
    # Convert results to MW and MVAR for printing
    P_g_MW = P_g / 1e6
    Q_total_MVAR = Q_total / 1e6
    Q_g_opt_MVAR = Q_g_opt / 1e6
    Q_opt_MVAR = Q_opt / 1e6

    # Print the step-by-step calculation with all numbers
    print("Step-by-Step Calculation of Optimal Reactive Power Injection:\n")
    print("1. System Parameters:")
    print(f"   Target Voltage (V_target): {V_target/1000} kV")
    print(f"   Line Impedance (Z_WF): {R_WF} + j{X_WF} Ohms")
    print(f"   Fault Resistance (R_F): {R_F} Ohms")
    print(f"   Generator Power Factor (PF_min): {PF_min} lagging\n")

    print("2. Required Power Calculation (to hold V_target):")
    print("   Total Active Power required from Bus-W (P_g):")
    print(f"   P_g = (V_target^2 * (R_WF + R_F)) / ((R_WF + R_F)^2 + X_WF^2)")
    print(f"   P_g = ({V_target}^2 * ({R_WF} + {R_F})) / (({R_WF} + {R_F})^2 + {X_WF}^2)")
    print(f"   P_g = {P_g_MW:.3f} MW\n")
    
    print("   Total Reactive Power required from Bus-W (Q_total):")
    print(f"   Q_total = (V_target^2 * X_WF) / ((R_WF + R_F)^2 + X_WF^2)")
    print(f"   Q_total = ({V_target}^2 * {X_WF}) / (({R_WF} + {R_F})^2 + {X_WF}^2)")
    print(f"   Q_total = {Q_total_MVAR:.3f} MVAR\n")
    
    print("3. Generator's Maximum Reactive Power (Q_g_max):")
    print(f"   Based on PF constraint (Q_g/P_g <= tan(acos({PF_min}))):")
    print(f"   Q_g_max = P_g * tan(acos(PF)) = {P_g_MW:.3f} MW * {tan_phi_max:.4f}")
    print(f"   Q_g_max = {Q_g_opt_MVAR:.3f} MVAR\n")

    print("4. Optimal Compensator Injection (Q_opt):")
    print("   To minimize Q_opt, the generator's reactive power output (Q_g) is maximized.")
    print(f"   Q_opt = Q_total - Q_g_max")
    print(f"   Q_opt = {Q_total_MVAR:.3f} MVAR - {Q_g_opt_MVAR:.3f} MVAR")
    print(f"   Q_opt = {Q_opt_MVAR:.3f} MVAR\n")
    
    # Check against compensator limit
    if Q_opt > Q_max:
        print(f"Warning: Required injection ({Q_opt_MVAR:.3f} MVAR) exceeds device max capacity ({Q_max/1e6} MVAR).")
    else:
        print(f"The optimal reactive power injection is within the device's capacity of {Q_max/1e6} MVAR.")

solve_reactive_power_optimization()