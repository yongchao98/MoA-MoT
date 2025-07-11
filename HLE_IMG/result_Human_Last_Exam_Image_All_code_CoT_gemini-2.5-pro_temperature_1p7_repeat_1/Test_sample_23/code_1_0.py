import math

def calculate_pump_work():
    """
    Calculates the work done by the pump based on the provided fluid dynamics problem.
    """
    # 1. Define all known parameters from the problem in SI units.
    rho = 997  # Density of water in kg/m^3
    p1 = 1.05e5  # Ambient pressure at tank surface (Point 1) in N/m^2
    pg = 1.33e5  # Absolute pressure at manometer (Point G) in N/m^2
    q = 2.86e-3  # Volume flow rate in m^3/s
    z1 = 2  # Height of water surface in tank (Point 1) in m
    zg = 3  # Height of manometer/mold inlet (Point G) in m
    r = 15.5 / 1000  # Pipe radius in m
    length = 14.9  # Pipe length in m
    f = 0.004  # Darcy friction factor
    k_shrinkage = 0.4  # Minor loss coefficient for tank exit (shrinkage)
    g = 9.81  # Acceleration due to gravity in m/s^2

    # Assumptions based on the problem description
    v1 = 0  # Velocity at the surface of the large tank is negligible

    print("--- Step 1: Initial Parameters ---")
    print(f"Density (ρ): {rho} kg/m^3")
    print(f"Pressure at Point 1 (P₁): {p1} N/m^2")
    print(f"Pressure at Point G (P_G): {pg} N/m^2")
    print(f"Volume Flow Rate (Q): {q} m^3/s")
    print(f"Elevation at Point 1 (z₁): {z1} m")
    print(f"Elevation at Point G (z_G): {zg} m")
    print(f"Pipe radius (r): {r} m")
    print(f"Pipe length (L): {length} m")
    print(f"Friction factor (f): {f}")
    print(f"Shrinkage loss coefficient (K_s): {k_shrinkage}")
    print(f"Gravity (g): {g} m/s^2\n")

    # 2. Calculate pipe properties and fluid velocity.
    d = 2 * r  # Pipe diameter
    area = math.pi * r**2  # Pipe cross-sectional area
    vg = q / area  # Velocity in the pipe (at Point G)
    vg_sq_over_2g = (vg**2) / (2 * g)  # Velocity head

    print("--- Step 2: Calculate Pipe & Flow Properties ---")
    print(f"Pipe Diameter (D): {d:.4f} m")
    print(f"Pipe Area (A): {area:.6f} m^2")
    print(f"Flow Velocity (v_G): {vg:.4f} m/s")
    print(f"Velocity Head (v_G²/2g): {vg_sq_over_2g:.4f} m\n")

    # 3. Calculate head loss from Point 1 to Point G.
    # This includes major loss (pipe friction) and minor loss (entrance).
    h_friction_1_to_g = (f * (length / d) + k_shrinkage) * vg_sq_over_2g

    print("--- Step 3: Calculate Head Loss (Point 1 to G) ---")
    print(f"Head Loss Formula: h_f = [f * (L/D) + K_s] * (v_G²/2g)")
    print(f"h_f = [{f} * ({length}/{d:.4f}) + {k_shrinkage}] * {vg_sq_over_2g:.4f}")
    print(f"h_f = {h_friction_1_to_g:.4f} m\n")

    # 4. Apply the energy equation from Point 1 to G to find pump head (h_pump).
    # h_pump = (P_G - P_1)/(ρ*g) + (v_G² - v_1²)/(2*g) + (z_G - z_1) + h_friction_1_to_g
    p_head_diff = (pg - p1) / (rho * g)
    v_head_diff = vg_sq_over_2g - (v1**2 / (2 * g))
    z_head_diff = zg - z1
    
    h_pump = p_head_diff + v_head_diff + z_head_diff + h_friction_1_to_g

    print("--- Step 4: Calculate Pump Head (h_pump) ---")
    print("Energy Equation: h_pump = (P_G - P₁)/(ρg) + (v_G² - v₁²)/(2g) + (z_G - z₁) + h_f")
    print(f"h_pump = ({pg} - {p1})/({rho}*{g}) + ({vg:.4f}² - {v1}²)/(2*{g}) + ({zg} - {z1}) + {h_friction_1_to_g:.4f}")
    print(f"h_pump = {p_head_diff:.4f} m + {v_head_diff:.4f} m + {z_head_diff} m + {h_friction_1_to_g:.4f} m")
    print(f"h_pump = {h_pump:.4f} m\n")
    
    # 5. Calculate the work of the pump (Power in Watts).
    # W_pump = ρ * g * Q * h_pump
    work_pump = rho * g * q * h_pump

    print("--- Step 5: Calculate Work of the Pump (Power) ---")
    print("Pump Power Formula: W_pump = ρ * g * Q * h_pump")
    print(f"W_pump = {rho} * {g} * {q} * {h_pump:.4f}")
    print(f"Final Pump Power (W_pump) = {work_pump:.2f} Watts\n")
    
    # Return final answer
    return work_pump

if __name__ == '__main__':
    final_answer = calculate_pump_work()
    # The final answer is requested at the end of the response.
    # print(f"<<<{final_answer:.2f}>>>")

calculate_pump_work()