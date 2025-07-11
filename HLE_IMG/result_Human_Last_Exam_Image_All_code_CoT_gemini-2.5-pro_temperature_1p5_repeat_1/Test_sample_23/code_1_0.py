import math

def calculate_pump_work():
    """
    Calculates the work of the pump based on the given fluid dynamics problem.
    """
    # 1. Define all constants and given values.
    rho = 997  # Density of water (kg/m^3)
    g = 9.81  # Acceleration due to gravity (m/s^2)
    Q = 2.86e-3  # Volume flow rate (m^3/s)
    z1 = 2  # Elevation of water surface in tank (m)
    z2 = 3  # Elevation of pipe exit (m)
    r_mm = 15.5  # Pipe radius (mm)
    L_D_ratio = 31  # L/D ratio for major loss calculation from table
    f = 0.004  # Darcy friction factor
    ef_shrinkage = 0.4  # Minor loss coefficient for contraction
    ef_expansion = 0.8  # Minor loss coefficient for expansion
    
    # Assumptions from the problem statement:
    # P1 = P2 = P_atm (Pressure at tank surface and exit are both ambient)
    # v1 = 0 (Velocity at the surface of a large tank is negligible)
    
    # 2. Convert units and calculate derived geometric properties.
    r_m = r_mm / 1000  # Convert radius from mm to m
    D = 2 * r_m  # Pipe diameter (m)
    A = math.pi * r_m**2  # Pipe cross-sectional area (m^2)
    
    # 3. Calculate fluid velocity in the pipe (v2).
    v = Q / A
    
    # 4. Calculate the kinetic head term.
    kinetic_head = (v**2) / (2 * g)
    
    # 5. Calculate head losses (h_L).
    # Major loss from pipe friction
    h_f = f * L_D_ratio * kinetic_head
    # Total minor loss coefficient
    K_L_total = ef_shrinkage + ef_expansion
    # Minor loss from fittings
    h_m = K_L_total * kinetic_head
    # Total head loss
    h_L = h_f + h_m
    
    # 6. Calculate the required pump head (h_pump) from the Bernoulli equation.
    # h_pump = (P2 - P1)/rho*g + (z2 - z1) + (v2^2 - v1^2)/2g + h_L
    # Since P1=P2 and v1=0, this simplifies to:
    h_pump = (z2 - z1) + kinetic_head + h_L
    
    # 7. Calculate the work of the pump (Power in Watts).
    W_pump = rho * g * Q * h_pump
    
    # 8. Print the calculation steps and the final answer.
    print("### Calculation of Pump Work (Power) ###")
    print("\nStep 1: Calculate fluid velocity and kinetic head.")
    print(f"Pipe radius: r = {r_mm} mm = {r_m:.4f} m")
    print(f"Pipe area: A = π * r² = {A:.6f} m²")
    print(f"Fluid velocity: v = Q / A = {Q} m³/s / {A:.6f} m² = {v:.3f} m/s")
    print(f"Kinetic head: v²/2g = ({v:.3f} m/s)² / (2 * {g} m/s²) = {kinetic_head:.3f} m")
    
    print("\nStep 2: Calculate head losses.")
    print(f"Major head loss: h_f = f * (L/D) * (v²/2g) = {f} * {L_D_ratio} * {kinetic_head:.3f} m = {h_f:.3f} m")
    print(f"Minor head loss: h_m = (K_shrinkage + K_expansion) * (v²/2g) = ({ef_shrinkage} + {ef_expansion}) * {kinetic_head:.3f} m = {h_m:.3f} m")
    print(f"Total head loss: h_L = h_f + h_m = {h_f:.3f} m + {h_m:.3f} m = {h_L:.3f} m")

    print("\nStep 3: Calculate the total head required from the pump.")
    print("h_pump = (z₂ - z₁) + (v²/2g) + h_L")
    print(f"h_pump = ({z2} m - {z1} m) + {kinetic_head:.3f} m + {h_L:.3f} m = {h_pump:.3f} m")

    print("\nStep 4: Calculate the final work of the pump (power).")
    print("W_pump = ρ * g * Q * h_pump")
    print(f"W_pump = {rho} kg/m³ * {g} m/s² * {Q} m³/s * {h_pump:.3f} m")
    print(f"W_pump = {W_pump:.2f} Watts")

calculate_pump_work()
<<<75.53>>>