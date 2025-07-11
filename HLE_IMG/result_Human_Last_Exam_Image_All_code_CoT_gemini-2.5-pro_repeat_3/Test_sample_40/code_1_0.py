import cmath
import math

def solve_hvac_optimization():
    """
    Solves the non-linear optimization problem for the HVAC transmission system.
    """
    # Step 1: Define given parameters in per-unit (p.u.)
    S_base = 100.0  # MVA
    V_base_B = 220.0 # kV
    
    # Transmission line impedance from Bus A to Bus B
    # Z_S = (0.02 + j0.10) Ohm is interpreted as p.u. as ohmic value is too small.
    Z_S = complex(0.02, 0.10) # p.u.
    
    # Fault impedance
    # Z_F = 0.15 Ohm is interpreted as a purely resistive p.u. value.
    Z_F = complex(0.15, 0) # p.u.
    
    # Voltages
    V_A = complex(1.0, 0)   # Voltage at Bus A (external grid), assumed as reference (1.0 p.u. at 0 degrees)
    V_B_target_mag = 1.0    # Target voltage magnitude at Bus B in p.u.

    # Maximum reactive power capacity of STATCOM
    Q_max_mvar = 50.0 # MVAR
    Q_max_pu = Q_max_mvar / S_base # p.u.

    print("--- System Parameters (in per-unit) ---")
    print(f"Base Power (S_base): {S_base} MVA")
    print(f"Line Impedance (Z_S): {Z_S:.4f} p.u.")
    print(f"Fault Impedance (Z_F): {Z_F:.4f} p.u.")
    print(f"Source Voltage (V_A): {V_A:.4f} p.u.")
    print(f"Target Voltage Magnitude at Bus B (|V_B|): {V_B_target_mag:.2f} p.u.")
    print("-" * 40)

    # Step 2: Formulate and solve the nodal equation
    # Admittances
    Y_S = 1 / Z_S
    Y_F = 1 / Z_F
    G_S, B_S = Y_S.real, Y_S.imag
    G_F, B_F = Y_F.real, Y_F.imag

    # The nodal equation at Bus B is: (V_B - V_A)*Y_S + V_B*Y_F + S_inj* / V_B* = 0
    # Let S_inj = jQ_comp. We want |V_B| = 1.0. Let V_B = 1.0 * exp(j*delta)
    # This leads to a real and imaginary part that must be zero.
    # Real part: G_S*cos(delta) + B_S*sin(delta) = G_S + G_F
    # Imaginary part: Q_comp = -B_S + B_S*cos(delta) - G_S*sin(delta)

    # Solve for delta from the real part equation: A*cos(d) + B*sin(d) = C
    A = G_S
    B = B_S
    C = G_S + G_F
    
    R_trig = math.sqrt(A**2 + B**2)
    alpha = cmath.phase(complex(A, B))
    
    # R*cos(delta - alpha) = C => cos(delta - alpha) = C / R
    cos_val = C / R_trig
    
    # We choose the stable solution for delta (smaller angle magnitude)
    delta_minus_alpha = math.acos(cos_val)
    delta = delta_minus_alpha + alpha # in radians

    delta_deg = math.degrees(delta)

    # Step 3: Calculate the required reactive power Q_opt (Q_comp in p.u.)
    Q_comp_pu = -B_S + B_S * math.cos(delta) - G_S * math.sin(delta)
    Q_opt_mvar = Q_comp_pu * S_base

    print("--- Optimization Results ---")
    print(f"Calculated voltage angle at Bus B (delta): {delta_deg:.2f} degrees")
    print(f"The minimum reactive power injection required is:")
    print(f"Q_opt = {Q_comp_pu:.4f} p.u. = {Q_opt_mvar:.2f} MVAR")
    if Q_opt_mvar > Q_max_mvar:
        print(f"Note: This required power ({Q_opt_mvar:.2f} MVAR) exceeds the STATCOM's maximum capacity of {Q_max_mvar} MVAR.")
    print("-" * 40)

    # Step 4: Calculate system power losses
    V_B = cmath.rect(V_B_target_mag, delta)
    
    # Current in the transmission line
    I_line = (V_A - V_B) / Z_S
    
    # Loss in the transmission line
    P_loss_line_pu = (abs(I_line)**2) * Z_S.real
    
    # Current through the fault
    I_fault = V_B / Z_F
    
    # Loss in the fault impedance
    P_loss_fault_pu = (abs(I_fault)**2) * Z_F.real
    
    # Total fundamental losses
    P_loss_base_pu = P_loss_line_pu + P_loss_fault_pu
    P_loss_base_mw = P_loss_base_pu * S_base
    
    # Step 5: Add harmonic losses
    harmonic_loss_factor = 0.04  # 4% increase
    P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_factor)
    P_loss_total_mw = P_loss_base_mw * (1 + harmonic_loss_factor)

    print("--- System Power Loss Calculation ---")
    print(f"Power loss in line (Z_S): {P_loss_line_pu * S_base:.2f} MW")
    print(f"Power loss in fault (Z_F): {P_loss_fault_pu * S_base:.2f} MW")
    print(f"Total fundamental loss: {P_loss_base_mw:.2f} MW")
    print(f"Total system loss including 4% harmonic effects = {P_loss_base_mw:.2f} MW * (1 + {harmonic_loss_factor}) = {P_loss_total_mw:.2f} MW")
    
    print("\n--- Final Answer ---")
    print(f"Optimal Reactive Power (Q_opt) = {Q_opt_mvar:.2f} MVAR")
    print(f"Total System Real Power Losses = {P_loss_total_mw:.2f} MW")
    
solve_hvac_optimization()