import math

def calculate_pump_work():
    """
    Calculates the work (power) of the pump based on the provided fluid dynamics parameters.
    """
    # Given parameters
    rho = 997  # density of water in kg/m^3
    P_ambient = 1.0e5  # ambient pressure in N/m^2 (Pa)
    P_m = 1.33e5  # manometer absolute pressure in N/m^2 (Pa)
    Q = 2.86e-3  # volume flow rate in m^3/s
    r_mm = 15.5  # pipe radius in mm
    g = 9.81  # acceleration due to gravity in m/s^2
    z1 = 2.0  # height of water in the tank above the pump in m
    ef_shrinkage = 0.4  # shrinkage (contraction) loss coefficient

    # --- Step 1: Calculate pipe properties ---
    r = r_mm / 1000.0  # convert radius to meters
    A = math.pi * r**2  # pipe cross-sectional area in m^2
    v_pipe = Q / A  # water velocity in the pipe in m/s

    # --- Step 2: Calculate pressure at pump inlet (P_in) using Bernoulli's equation ---
    # from tank surface (1) to pump inlet (in).
    # P_in/rho = P_ambient/rho + g*z1 - v_pipe^2/2 - w_loss_inlet
    # w_loss_inlet = ef_shrinkage * (v_pipe**2 / 2)
    # P_in = P_ambient + rho * (g*z1 - (1 + ef_shrinkage) * v_pipe**2 / 2)
    
    ke_term = v_pipe**2 / 2 # Kinetic energy per unit mass
    pressure_head_at_inlet = P_ambient + rho * (g * z1 - (1 + ef_shrinkage) * ke_term)
    P_in = pressure_head_at_inlet
    
    # --- Step 3: Calculate work per unit mass (w_p) across the pump ---
    # w_p = (P_out - P_in) / rho, where P_out = P_m
    w_p = (P_m - P_in) / rho

    # --- Step 4: Calculate the total power of the pump ---
    mass_flow_rate = rho * Q
    power = w_p * mass_flow_rate

    # --- Output the results ---
    print("--- Calculation Steps ---")
    print(f"1. Pipe velocity (v) = {v_pipe:.3f} m/s")
    print(f"2. Pressure at pump inlet (P_in) = {P_in:.2f} N/m^2")
    print(f"3. Work per unit mass (w_p) = {w_p:.3f} J/kg")
    
    print("\n--- Final Calculation for Pump Power (Work) ---")
    print(f"Equation: Power = (w_p) * (mass_flow_rate)")
    print(f"Power = ({w_p:.3f}) * ({mass_flow_rate:.3f})")
    print(f"Work of the pump = {power:.1f} Watts")
    
    return power

# Run the calculation and store the result
final_work = calculate_pump_work()
print(f"\n<<<67.1>>>")
