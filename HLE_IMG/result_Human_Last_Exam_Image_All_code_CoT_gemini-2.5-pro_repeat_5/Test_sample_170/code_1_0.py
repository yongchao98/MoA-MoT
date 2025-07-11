import cmath
import math
import numpy as np

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """
    # 1. System Parameters and Per-Unit Conversion
    S_base_MVA = 100.0  # Base Power in MVA
    V_nom_B_kV = 220.0   # Nominal Voltage at Bus B in kV
    
    # Fault condition: Voltage at Bus B drops to 85%
    V_fault_pu = 0.85
    
    # Target voltage at Bus B is nominal voltage
    V_target_pu = 1.0
    
    # System impedance (assumed to be Thevenin equivalent in pu)
    # Z_th = Z_S = (0.02 + j0.10) pu
    Z_th_real = 0.02
    Z_th_imag = 0.10
    Z_th = complex(Z_th_real, Z_th_imag)
    
    # Harmonic loss increase factor
    harmonic_loss_factor = 0.04
    
    # Maximum STATCOM capacity (for context, not used in calculation as we find the required Q)
    Q_max_MVAR = 50.0
    Q_max_pu = Q_max_MVAR / S_base_MVA

    # 2. Formulate and Solve the Quadratic Equation for Q_C
    # The relationship between voltage and reactive power injection Q_C is:
    # |Z_th|^2 * Q_C^2 - (2 * |V_B|^2 * X_th) * Q_C + |V_B|^2 * (|V_B|^2 - |V_th|^2) = 0
    # where |V_th| is the initial voltage (V_fault_pu) and |V_B| is the final voltage (V_target_pu)
    
    V_th_pu = V_fault_pu # Thevenin voltage is the voltage at the bus before compensation
    
    # Coefficients of the quadratic equation a*x^2 + b*x + c = 0 for Q_C
    a = abs(Z_th)**2
    b = -2 * V_target_pu**2 * Z_th.imag
    c = V_target_pu**2 * (V_target_pu**2 - V_th_pu**2)
    
    # Solve the quadratic equation
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("No real solutions for Q_C exist. Voltage restoration to 1.0 pu is not possible.")
        return

    sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
    sol2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    # 3. Determine the minimum required reactive power (Q_opt)
    # We are looking for the minimum positive injection of reactive power.
    solutions = [s for s in [sol1, sol2] if s > 0]
    
    if not solutions:
        print("No positive reactive power injection solution found.")
        return
        
    Q_opt_pu = min(solutions)
    Q_opt_MVAR = Q_opt_pu * S_base_MVA

    # 4. Calculate System Real Power Losses
    # Loss is due to current flowing through the Thevenin resistance
    # I_C^2 = Q_opt_pu^2 / V_target_pu^2
    I_C_sq = Q_opt_pu**2 / V_target_pu**2
    
    # Base power loss in pu
    P_loss_pu = I_C_sq * Z_th.real
    
    # Total loss including 4% harmonic increase
    P_loss_total_pu = P_loss_pu * (1 + harmonic_loss_factor)
    P_loss_total_MW = P_loss_total_pu * S_base_MVA

    # 5. Print the results
    print("--- Optimization Results ---")
    print(f"The minimum reactive power injection required from the STATCOM is:")
    print(f"Q_opt = {Q_opt_pu:.3f} pu")
    print(f"This is equivalent to: {Q_opt_MVAR:.2f} MVAR\n")
    
    print("--- System Power Loss Calculation ---")
    print(f"The total system real power losses after compensation are:")
    print(f"P_loss = {P_loss_total_pu:.4f} pu")
    print(f"This is equivalent to: {P_loss_total_MW:.2f} MW\n")

    print("--- Final Equation for Q_opt ---")
    print("The optimal reactive power Q_opt is the minimum positive root of the quadratic equation:")
    print(f"({a:.4f}) * Q_C^2 + ({b:.2f}) * Q_C + ({c:.4f}) = 0")


solve_hvac_optimization()