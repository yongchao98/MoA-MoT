import cmath
import math

def solve_wind_energy_optimization():
    """
    Solves the nonlinear optimization problem for the wind energy system.
    """
    # 1. Define constants in per-unit (p.u.)
    V_W_target = 0.575  # Target voltage at Bus-W (p.u.)
    R_WF = 0.01         # Transmission line resistance (p.u.)
    X_WF = 0.05         # Transmission line reactance (p.u.)
    R_F = 0.1           # Fault resistance (p.u.)
    S_base = 100        # Base power (MVA)
    Q_max_MVAR = 10     # Max reactive power of compensator (MVAR)
    Q_max_pu = Q_max_MVAR / S_base # Max reactive power (p.u.)
    
    # Harmonic loss factor (1 + 2% + 4% = 1.06)
    loss_factor = 1.06
    
    # Power factor constraint
    PF = 0.95
    
    # 2. Calculate helper constants
    # Relationship between Q and P from the power factor: Q = k * P
    k = math.tan(math.acos(PF))
    
    # 3. Set up and solve the quadratic equation for P_W
    # The system model leads to a quadratic equation for P_W of the form a*P_W^2 + b*P_W + c = 0
    # Based on the power balance and voltage drop equations:
    # P_W = P_fault + P_loss
    # P_fault = |V_F|^2 / R_F
    # |V_F| approx V_W - (P_W*R_WF + Q_W*X_WF)/V_W
    # P_loss = loss_factor * |I_W|^2 * R_WF
    
    # Coefficients of the quadratic equation a*x^2 + b*x + c = 0 where x = P_W
    R_plus_kX = R_WF + k * X_WF
    
    a = (R_plus_kX)**2 + loss_factor * R_F * R_WF * (1 + k**2)
    b = -2 * V_W_target**2 * R_plus_kX - R_F * V_W_target**2
    c = V_W_target**4
    
    # Solve the quadratic equation for P_W
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("No real solution for P_W exists. The system cannot be stabilized under these conditions.")
        return

    P_W1 = (-b + math.sqrt(discriminant)) / (2 * a)
    P_W2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    # Choose the physically meaningful (smaller) solution for power
    P_W = min(P_W1, P_W2)
    if P_W < 0:
         P_W = max(P_W1, P_W2)


    # 4. Calculate the required reactive power injection Q_comp
    # From reactive power balance: Q_W = Q_comp + Q_fault_side
    # Q_comp = Q_W - Q_loss
    # Q_W = k * P_W
    # Q_loss = loss_factor * |I_W|^2 * X_WF
    
    I_W_squared = (P_W**2 * (1 + k**2)) / V_W_target**2
    Q_loss = loss_factor * I_W_squared * X_WF
    Q_W_gen = k * P_W
    
    Q_opt_pu = Q_loss - Q_W_gen
    
    # Convert result to MVAR
    Q_opt_MVAR = Q_opt_pu * S_base
    
    # 5. Print the results and the final equation
    print("--- System Parameters and Assumptions ---")
    print(f"Target Voltage at Bus-W (V_W): {V_W_target} p.u.")
    print(f"Line Impedance (Z_WF): {R_WF} + j{X_WF} p.u.")
    print(f"Fault Resistance (R_F): {R_F} p.u.")
    print(f"Power Factor (PF): {PF} lagging")
    print(f"Harmonic Loss Factor: {loss_factor}")
    print("\n--- Calculation Steps ---")
    print(f"1. Solved for active power from wind farm (P_W): {P_W:.4f} p.u.")
    
    print("\n2. The optimal reactive power injection (Q_opt) is calculated using the reactive power balance equation:")
    print("   Q_opt = (Total Reactive Power Loss) - (Reactive Power from Wind Farm)")
    
    # Print the final equation with numerical values
    equation_str = (
        f"   Q_opt = ( {loss_factor:.2f} * (({P_W:.4f}^2 * (1 + {k:.4f}^2)) / {V_W_target:.3f}^2) * {X_WF} ) - ( {k:.4f} * {P_W:.4f} )"
    )
    print(equation_str)
    
    # Breakdown of the equation
    print(f"   Q_opt = ( {loss_factor:.2f} * ({I_W_squared:.4f}) * {X_WF} ) - ( {Q_W_gen:.4f} )")
    print(f"   Q_opt = ( {Q_loss:.4f} ) - ( {Q_W_gen:.4f} )")
    
    print("\n--- Final Result ---")
    print(f"The optimal reactive power injection required is {Q_opt_pu:.4f} p.u.")
    print(f"This corresponds to {Q_opt_MVAR:.2f} MVAR.")
    
    if Q_opt_pu > Q_max_pu:
        print(f"\nWarning: The required reactive power ({Q_opt_MVAR:.2f} MVAR) exceeds the device's maximum capacity of {Q_max_MVAR} MVAR.")
        print("This indicates that the voltage target cannot be met without relaxing other constraints (e.g., power factor).")

solve_wind_energy_optimization()