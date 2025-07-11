import math

def solve_pump_work():
    """
    Calculates the work of the pump based on the given system parameters.
    """

    # --- Given values and constants ---
    # Density of water in kg/m^3
    rho = 997
    # Volume flow rate in m^3/s
    Q = 2.86e-3
    # Height of water in tank above the pump's level in meters
    z_tank = 2.0
    # Height of pipe exit above the pump's level in meters
    z_exit = 3.0
    # Total length of the pipe in meters
    L_total = 14.9
    # Pipe radius in meters
    r = 15.5 / 1000
    # Darcy friction factor (dimensionless)
    f = 0.004
    # Minor loss coefficient for entrance shrinkage (dimensionless)
    ef_shrink = 0.4
    # Minor loss coefficient for exit expansion (dimensionless)
    ef_expand = 0.8
    # Acceleration due to gravity in m/s^2
    g = 9.81
    
    # Assumption: The provided ambient pressure "105 N/m^2" is a typo for 1.0x10^5 N/m^2.
    P_atm = 1.0e5
    # Assumption: The manometer pressure "1.33 x 105 N/m^2" is a typo for 1.33x10^5 N/m^2.
    P_p = 1.33e5

    # --- Step 1: Calculate basic pipe and flow properties ---
    D = 2 * r  # Pipe diameter, m
    A = math.pi * r**2  # Pipe cross-sectional area, m^2
    v_pipe = Q / A  # Fluid velocity in the pipe, m/s
    v_head = v_pipe**2 / (2 * g)  # Velocity head, m

    # --- Step 2: Analyze discharge side (manometer to exit) to find head loss ---
    # Energy Equation: (P_p/ρg) + v_p²/2g + z_p = (P_atm/ρg) + v_exit²/2g + z_exit + h_L_discharge
    # Assumptions: z_p=0 (datum), v_exit≈0.
    h_L_discharge = (P_p - P_atm) / (rho * g) + v_head - z_exit

    # --- Step 3: Use discharge head loss to find discharge pipe length ---
    # h_L_discharge = (f * L_discharge/D + ef_expand) * v_head
    L_discharge = ((h_L_discharge / v_head) - ef_expand) * D / f

    # --- Step 4: Calculate suction pipe length and suction side head loss ---
    L_suction = L_total - L_discharge
    h_L_suction = (f * L_suction / D + ef_shrink) * v_head

    # --- Step 5: Calculate total head loss and the required pump head (Hp) ---
    h_L_total = h_L_suction + h_L_discharge
    # Overall energy balance (tank to exit): H_p = (z_exit - z_tank) + h_L_total
    Hp = (z_exit - z_tank) + h_L_total

    # --- Step 6: Calculate the work of the pump (Power) in Watts ---
    W_pump = Hp * rho * g * Q

    # --- Step 7: Print the final calculation and result ---
    print("Final Calculation of Pump Work (Power):")
    print("Work = Pump_Head * Density * Gravity * Flow_Rate")
    print(f"Work = {Hp:.2f} m * {rho} kg/m^3 * {g} m/s^2 * {Q:.4f} m^3/s")
    print(f"\nThe work of the pump is {W_pump:.1f} Watts.")
    
    return W_pump

# Execute the function
solve_pump_work()