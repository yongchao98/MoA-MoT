import math

def calculate_pump_work():
    """
    Calculates the work of the pump based on the provided system parameters.
    """
    # Given values
    rho = 997  # Density of water in kg/m^3
    P1 = 105  # Ambient pressure at the source in N/m^2
    P_m = 1.33e5  # Absolute pressure at the manometer in N/m^2
    Q = 2.86e-3  # Volume flow rate in m^3/s
    z1 = 2  # Height of water in the tank in meters
    z_m = 0 # Height of the manometer (datum)
    r = 15.5 / 1000  # Pipe radius in meters
    L = 14.9  # Total pipe length in meters
    f = 0.004  # Darcy friction factor
    K_ent = 0.4  # Loss coefficient for entrance (shrinkage)
    LD_fitting = 31  # L/D for one bend
    num_bends = 2 # Number of bends
    g = 9.81  # Acceleration due to gravity in m/s^2

    # Step 1: Calculate pipe diameter, area, velocity, and velocity head
    D = 2 * r
    A = math.pi * r**2
    v = Q / A
    v_sq_2g = v**2 / (2 * g)

    # Step 2: Calculate head loss from tank to manometer (H_L_1_m)
    # Major loss (pipe friction)
    K_pipe = f * (L / D)
    # Minor losses (entrance + bends)
    # The loss for bends is calculated using the equivalent length method: K = f * (L/D)_eq
    K_bends = num_bends * f * LD_fitting
    # Total loss coefficient
    K_total = K_ent + K_bends + K_pipe
    # Total head loss
    H_L_1_m = K_total * v_sq_2g

    # Step 3: Apply the energy equation to find the pump head (H_p)
    # P1/(rho*g) + z1 + H_p = P_m/(rho*g) + z_m + v^2/(2*g) + H_L_1_m
    # H_p = (P_m - P1)/(rho*g) + (z_m - z1) + v^2/(2*g) + H_L_1_m
    term_pressure_head = (P_m - P1) / (rho * g)
    term_elevation_head = z_m - z1
    
    H_p = term_pressure_head + term_elevation_head + v_sq_2g + H_L_1_m

    # Step 4: Calculate the specific work of the pump (w_p = g * H_p)
    # w_p = g * H_p = (P_m - P1)/rho + g*(z_m - z1) + v^2/2 + g*H_L_1_m
    # A simpler form is w_p = (P_m - P1)/rho + g*(z_m - z1) + (1 + K_total) * v^2/2
    w_p = g * H_p

    # Print the results and the final equation
    print("--- Calculation Steps ---")
    print(f"Pipe diameter D = {D:.4f} m")
    print(f"Pipe area A = {A:.6f} m^2")
    print(f"Water velocity v = {v:.4f} m/s")
    print(f"Kinetic head v^2/(2g) = {v_sq_2g:.4f} m")
    print("\n--- Head Loss Calculation (from tank to manometer) ---")
    print(f"Pipe friction loss coefficient (f*L/D) = {K_pipe:.4f}")
    print(f"Bends loss coefficient (2*f*L/D_fitting) = {K_bends:.4f}")
    print(f"Entrance loss coefficient = {K_ent:.4f}")
    print(f"Total loss coefficient K_total = {K_total:.4f}")
    print(f"Total head loss H_L = K_total * v^2/(2g) = {H_L_1_m:.4f} m")
    print("\n--- Pump Head Calculation ---")
    print(f"Pressure head term = {term_pressure_head:.4f} m")
    print(f"Elevation head term = {term_elevation_head:.4f} m")
    print(f"Pump Head H_p = {H_p:.4f} m")
    print("\n--- Final Pump Work Calculation ---")
    print(f"Specific work of the pump w_p = g * H_p = {w_p:.2f} J/kg")
    
    print("\n--- Final Equation with numerical values ---")
    # Equation for specific work: w_p = g * H_p
    # w_p = g * [ (P_m - P1)/(rho*g) + (z_m - z1) + v^2/2g + H_L_1_m ]
    print(f"Work = {g:.2f} * [ ({P_m:.0f} - {P1:.0f})/({rho:.0f}*{g:.2f}) + ({z_m:.0f} - {z1:.0f}) + {v_sq_2g:.4f} + {H_L_1_m:.4f} ]")
    print(f"Work = {g:.2f} * [ {term_pressure_head:.4f} + {term_elevation_head:.1f} + {v_sq_2g:.4f} + {H_L_1_m:.4f} ]")
    print(f"Work = {g:.2f} * [ {H_p:.4f} ] = {w_p:.2f} J/kg")


calculate_pump_work()
<<<139.31>>>