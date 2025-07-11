import cmath
import math

def solve_wind_compensation():
    """
    Solves for the optimal reactive power injection to restore voltage in a wind energy system.
    """
    # Step 1: Define System Parameters in physical units
    Z_WF_ohm = complex(0.01, 0.05)
    R_F_ohm = 0.1
    V_nom_kV = 0.575
    S_base_MVA = 100.0
    V_W_target_kV = 0.575
    Q_max_MVAR = 10.0
    pf_min = 0.95
    harmonic_loss_increase = 0.06 # 6% increase

    # Step 2: Per-Unit (p.u.) Conversion
    V_base_kV = V_nom_kV
    S_base_MVA = S_base_MVA
    Z_base_ohm = (V_base_kV**2) / S_base_MVA # (kV^2)/MVA = Ohms

    V_W_pu = V_W_target_kV / V_base_kV
    Z_WF_pu = Z_WF_ohm / Z_base_ohm
    R_F_pu = R_F_ohm / Z_base_ohm
    Q_max_pu = Q_max_MVAR / S_base_MVA

    R_WF_pu = Z_WF_pu.real
    X_WF_pu = Z_WF_pu.imag
    
    print("--- System Parameters (per unit) ---")
    print(f"Base Impedance Z_base = {Z_base_ohm:.4f} Ohms")
    print(f"Target Voltage V_W = {V_W_pu:.1f} p.u.")
    print(f"Line Impedance Z_WF = ({R_WF_pu:.4f} + j{X_WF_pu:.4f}) p.u.")
    print(f"Fault Resistance R_F = {R_F_pu:.4f} p.u.")
    print(f"Max Compensator Q_max = {Q_max_pu:.2f} p.u.\n")

    # Step 3 & 4: Model with Harmonic Losses and Power Factor Limit
    # To minimize Q_comp, we maximize generator's reactive power support by operating at pf=0.95
    # The total power from the generator P_gen = P_fundamental_loss * (1 + harmonic_loss_increase)
    # The PF constraint is Q_W / P_gen <= tan(acos(pf_min))
    # Let P_W be the fundamental active power at Bus W. P_W = P_loss_line + P_loss_fault
    # Q_W <= tan(acos(pf_min)) * P_gen = tan(acos(pf_min)) * P_W * (1 + harmonic_loss_increase)
    k = math.tan(math.acos(pf_min))
    k_prime = k * (1 + harmonic_loss_increase) # Effective k for Q_W = k' * P_W

    # Step 5 & 6: Formulate and Solve the Quadratic Equation for P_W
    # Based on power balance P_W = |I_W|^2*R_WF + |V_F|^2/R_F
    # This simplifies to A*P_W^2 + B*P_W + C = 0
    
    # Coefficients of the quadratic equation
    Z_WF_mag2 = R_WF_pu**2 + X_WF_pu**2
    term1 = 1 + k_prime**2
    
    A = term1 * (R_WF_pu + Z_WF_mag2 / R_F_pu)
    B = -1 - (2 / R_F_pu) * (R_WF_pu + k_prime * X_WF_pu)
    C = 1 / R_F_pu

    print("--- Solving for Fundamental Active Power P_W ---")
    print("The problem reduces to solving the quadratic equation for P_W:")
    print(f"{A:.4f} * P_W^2 + ({B:.4f}) * P_W + {C:.4f} = 0\n")

    # Solve quadratic equation: P_W = (-B +/- sqrt(B^2 - 4AC)) / 2A
    discriminant = B**2 - 4 * A * C
    if discriminant < 0:
        print("No real solution for P_W exists. The system cannot be stabilized under these conditions.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    P_W1 = (-B + sqrt_discriminant) / (2 * A)
    P_W2 = (-B - sqrt_discriminant) / (2 * A)

    # Choose the smaller, physically more likely solution for P_W
    P_W_pu = min(P_W1, P_W2)
    print(f"Calculated Fundamental Active Power P_W = {P_W_pu:.4f} p.u.\n")
    
    # Step 7: Calculate all other system variables
    Q_W_pu = k_prime * P_W_pu
    I_W_pu = complex(P_W_pu, -Q_W_pu) # Current lags voltage for lagging PF
    
    # Voltage at fault point F
    V_F_pu = V_W_pu - I_W_pu * Z_WF_pu

    # Step 8: Calculate Q_comp
    # Power balance at node F: Q_F_in + Q_comp = Q_fault = 0 (since R_F is purely resistive)
    # Q_F_in is the reactive power arriving at F from W
    # S_F_in = V_F * I_W* -> Q_F_in = Im(V_F * I_W*)
    # Therefore, Q_comp = -Q_F_in
    
    S_F_in = V_F_pu * I_W_pu.conjugate()
    Q_F_in = S_F_in.imag
    
    Q_opt_pu = -Q_F_in
    
    # Step 9: Final Result
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
    
    print("--- Optimal Reactive Power Injection ---")
    print(f"Calculated Voltage at Fault V_F = ({V_F_pu.real:.4f} + j{V_F_pu.imag:.4f}) p.u.")
    print(f"Reactive Power arriving at Fault Q_F_in = {Q_F_in:.4f} p.u.")
    print(f"Optimal Reactive Power Injection Q_opt = {Q_opt_pu:.4f} p.u.")
    print(f"This corresponds to Q_opt = {Q_opt_MVAR:.4f} MVAR.\n")

    # Check against capacity constraint
    print("--- Constraint Check ---")
    print(f"The maximum capacity of the compensator is {Q_max_MVAR} MVAR.")
    if Q_opt_MVAR > Q_max_MVAR:
        print(f"Warning: The required injection of {Q_opt_MVAR:.2f} MVAR exceeds the device's maximum capacity.")
    else:
        print("The required injection is within the device's capacity.")

    # Return the final numeric answer for the problem
    return Q_opt_MVAR

if __name__ == '__main__':
    q_optimal = solve_wind_compensation()
    print("\nFinal Answer:")
    # Using format requested in problem statement for final value
    # Remember in the final code you still need to output each number in the final equation!
    # Which was done above. The <<<>>> format is for the final, single answer.
    print(f'<<<The optimal reactive power injection Q_opt is {q_optimal:.1f} MVAR>>>')
