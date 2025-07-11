import math

def solve_reactive_power_optimization():
    """
    Solves for the optimal reactive power injection to restore voltage in a wind energy system.
    """

    # 1. Define System Parameters in Per Unit (p.u.) and Physical Values
    S_base_MVA = 100.0  # MVA
    V_base_transmission_kV = 400.0  # kV
    V_W_target_kV = 0.575 # This is stated for the generator bus, but target for Bus-W
                           # Assuming this is a typo and means 0.575 p.u.
    
    # Transmission line impedance (from Bus-W to Fault F)
    R_WF_ohm = 0.01
    X_WF_ohm = 0.05
    
    # Fault resistance
    R_F_ohm = 0.1
    
    # Maximum reactive power from the compensator
    Q_max_MVAR = 10.0
    
    # Wind farm operational parameters
    PF_min_lagging = 0.95
    # Assume wind farm active power output is 90 MW (a reasonable assumption)
    P_farm_MW = 90.0

    # Harmonic loss factor
    harmonic_loss_increase = 0.06 # 6% increase (2% + 4%)
    loss_factor = 1.0 + harmonic_loss_increase
    
    # --- Convert to Per Unit ---
    V_base_bus_W = V_base_transmission_kV 
    Z_base_ohm = V_base_bus_W**2 / S_base_MVA

    # Assuming given impedances are already in p.u. as ohmic values lead to physically inconsistent results.
    R_WF_pu = R_WF_ohm 
    X_WF_pu = X_WF_ohm
    R_F_pu = R_F_ohm

    V_W_target_pu = 0.575 # Target voltage at Bus-W in p.u.
    V_G_pu = 1.0 # Grid voltage (infinite bus) in p.u.
    
    P_W_pu = P_farm_MW / S_base_MVA

    # Calculate reactive power support from the wind farm (Q_W)
    # The farm operates at its maximum reactive power capability to support voltage
    phi = math.acos(PF_min_lagging)
    Q_W_pu = P_W_pu * math.tan(phi)
    
    # 2. Model the System and Incorporate Losses
    # We model the system as a two-bus model (Grid and Bus-W) connected by an equivalent impedance.
    # The harmonic losses are modeled by increasing the resistance of the line segment Z_WF.
    R_WF_effective_pu = R_WF_pu * loss_factor
    
    # The equivalent impedance Z_eq is modeled by summing the line impedance to the fault 
    # and the fault resistance. This is a simplification to make the problem solvable.
    R_eq_pu = R_WF_effective_pu + R_F_pu
    X_eq_pu = X_WF_pu
    Z_eq_mag = math.sqrt(R_eq_pu**2 + X_eq_pu**2)
    Z_eq_angle_rad = math.atan2(X_eq_pu, R_eq_pu)

    # 3. Solve Power Flow Equations
    # The optimization goal min(Q_c) with V_W=constant is achieved by maximizing P_W.
    # We use P_W = 0.9 p.u. and solve for the required Q_total = Q_W + Q_c.
    
    # Power flow equation for Active Power (P_W):
    # P_W = (|V_G||V_W|/|Z_eq|) * cos(theta_z + delta) - (|V_W|^2/|Z_eq|) * cos(theta_z)
    
    # --- Solve for voltage angle delta ---
    # Rearranging the equation for cos(theta_z + delta):
    # cos(theta_z + delta) = (P_W + (|V_W|^2/|Z_eq|) * cos(theta_z)) / (|V_G||V_W|/|Z_eq|)
    
    A = (V_G_pu * V_W_target_pu) / Z_eq_mag
    B = (V_W_target_pu**2) / Z_eq_mag
    
    cos_theta_z = math.cos(Z_eq_angle_rad)
    
    # Equation for cos(theta_z + delta)
    cos_angle_arg = (P_W_pu + B * cos_theta_z) / A
    
    print("--- Solving for Optimal Reactive Power Injection ---")
    print("\nStep 1: System Parameters (p.u.)")
    print(f"  P_W (Wind Farm Active Power) = {P_W_pu:.3f} p.u.")
    print(f"  Q_W (Wind Farm Reactive Power) = {Q_W_pu:.3f} p.u.")
    print(f"  V_W (Target Voltage at Bus-W) = {V_W_target_pu:.3f} p.u.")
    print(f"  V_G (Grid Voltage) = {V_G_pu:.3f} p.u.")
    print(f"  Z_eq (Equivalent Impedance) = {R_eq_pu:.4f} + j{X_eq_pu:.4f} p.u.")
    print(f"  |Z_eq| = {Z_eq_mag:.4f} p.u., Angle(Z_eq) = {math.degrees(Z_eq_angle_rad):.2f} degrees")
    print(f"  Harmonic loss factor applied to R_WF: {loss_factor:.2f}")


    print("\nStep 2: Solving Active Power Equation for Bus-W angle (delta)")
    print(f"  P_W = (V_G*V_W/|Z_eq|) * cos(theta_z + delta) - (V_W^2/|Z_eq|) * cos(theta_z)")
    print(f"  {P_W_pu:.3f} = ({A:.3f}) * cos({Z_eq_angle_rad:.3f} + delta) - ({B:.3f}) * cos({Z_eq_angle_rad:.3f})")
    
    # Handle case where cos_angle_arg might be out of [-1, 1] due to unrealistic parameters
    if abs(cos_angle_arg) > 1:
        print("\nError: The specified power transfer is impossible under the given voltage constraints.")
        print("The system parameters lead to a situation with no real solution for the voltage angle.")
        return

    angle_rad = math.acos(cos_angle_arg)
    
    # Two possible solutions for angle, we take the one that corresponds to stable operation
    delta_rad = angle_rad - Z_eq_angle_rad

    print(f"  cos({Z_eq_angle_rad:.3f} + delta) = {cos_angle_arg:.4f}")
    print(f"  --> {Z_eq_angle_rad:.3f} + delta = {angle_rad:.4f} rad")
    print(f"  --> delta = {delta_rad:.4f} rad ({math.degrees(delta_rad):.2f} degrees)")

    # --- Solve for total reactive power Q_tot ---
    # Q_tot = (|V_G||V_W|/|Z_eq|) * sin(theta_z + delta) - (|V_W|^2/|Z_eq|) * sin(theta_z)
    sin_theta_z = math.sin(Z_eq_angle_rad)
    sin_angle_arg = math.sin(angle_rad)
    
    Q_tot_pu = A * sin_angle_arg - B * sin_theta_z

    print("\nStep 3: Solving Reactive Power Equation for Q_total")
    print(f"  Q_total = (V_G*V_W/|Z_eq|) * sin(theta_z + delta) - (V_W^2/|Z_eq|) * sin(theta_z)")
    print(f"  Q_total = ({A:.3f}) * sin({angle_rad:.4f}) - ({B:.3f}) * sin({Z_eq_angle_rad:.3f})")
    print(f"  Q_total = {Q_tot_pu:.4f} p.u.")
    
    # --- Calculate the compensator's reactive power Q_c ---
    Q_c_pu = Q_tot_pu - Q_W_pu
    Q_c_MVAR = Q_c_pu * S_base_MVA
    
    print("\nStep 4: Final Calculation of Q_opt")
    print(f"  Q_opt = Q_total - Q_W")
    print(f"  Q_opt = {Q_tot_pu:.4f} p.u. - {Q_W_pu:.4f} p.u. = {Q_c_pu:.4f} p.u.")
    print(f"\n--- Result ---")
    print(f"The optimal reactive power injection Q_opt required to restore the voltage at Bus-W to {V_W_target_pu} p.u. is:")
    print(f"Q_opt = {Q_c_MVAR:.1f} MVAR")
    
    if Q_c_MVAR > Q_max_MVAR:
        print(f"\nWarning: The required reactive power ({Q_c_MVAR:.1f} MVAR) exceeds the device's maximum capacity of {Q_max_MVAR:.1f} MVAR.")
        print("The voltage restoration target cannot be fully achieved with the available compensation.")
        
    return Q_c_MVAR

# Execute the function to find the solution.
q_optimal_mvar = solve_reactive_power_optimization()
# Print the final numerical answer in the required format.
# Note: The calculation shows the target is physically unreachable with the given device, 
# but the question asks for the optimal injection *to restore the voltage*, which we have calculated.
# The optimal feasible injection would be Q_max, but that would not meet the voltage target.
# Thus we provide the calculated value.
final_answer = round(q_optimal_mvar,1)
print(f"\n<<< {final_answer} >>>")