import math

def calculate_pump_work():
    """
    Calculates the work of the pump based on the provided fluid dynamics problem.
    """
    # --- Given values from the problem description and table ---
    rho = 997  # Density of water in kg/m^3
    g = 9.81  # Acceleration of gravity in m/s^2
    Q = 2.86e-3  # Volume flow rate in m^3/s
    
    # Heights relative to the pump
    z1 = 2.0  # Height of water in the tank in meters
    z2 = 3.0  # Vertical length of the pipe to the exit in meters
    
    # Pipe and friction parameters
    r_mm = 15.5  # Pipe radius in mm
    L = 14.9  # Total pipe length in m
    f = 0.004  # Darcy friction factor
    LD_fittings = 31  # L/D ratio for fittings
    K_entrance = 0.4  # Loss coefficient for shrinkage (entrance)
    K_exit = 0.8  # Loss coefficient for expansion (exit)
    
    # --- Step 1: Calculate geometric and flow properties ---
    r = r_mm / 1000  # Convert radius to meters
    D = 2 * r  # Pipe diameter in meters
    A = math.pi * r**2  # Pipe cross-sectional area in m^2
    v = Q / A  # Water velocity in the pipe in m/s
    
    v_sq_over_2 = v**2 / 2  # Kinetic energy term per unit mass in J/kg

    print("--- Intermediate Calculations ---")
    print(f"Pipe Diameter (D): {D:.4f} m")
    print(f"Pipe Area (A): {A:.6f} m^2")
    print(f"Flow Velocity (v): {v:.4f} m/s")
    print(f"Kinetic Energy per mass (v^2/2): {v_sq_over_2:.4f} J/kg\n")

    # --- Step 2: Calculate total friction loss (W_friction_total) ---
    # Major loss coefficient
    K_major = f * (L / D)
    
    # Minor loss coefficient for fittings
    K_fittings = f * LD_fittings

    # Total equivalent loss coefficient
    K_total = K_major + K_entrance + K_fittings + K_exit
    
    # Total friction loss per unit mass
    W_friction_total = K_total * v_sq_over_2

    print("--- Friction Loss Calculation ---")
    print(f"Major Loss Coefficient (f*L/D): {K_major:.4f}")
    print(f"Minor Loss Coefficient (Entrance): {K_entrance:.4f}")
    print(f"Minor Loss Coefficient (Fittings): {K_fittings:.4f}")
    print(f"Minor Loss Coefficient (Exit): {K_exit:.4f}")
    print(f"Total Loss Coefficient (K_total): {K_total:.4f}")
    print(f"Total Friction Loss per mass (W_friction): {W_friction_total:.4f} J/kg\n")
    
    # --- Step 3: Calculate potential energy change ---
    # Potential energy change per unit mass
    W_potential = g * (z2 - z1)

    print("--- Potential Energy Calculation ---")
    print(f"Potential Energy Change per mass (g*(z2-z1)): {W_potential:.4f} J/kg\n")

    # --- Step 4: Final calculation for pump work ---
    # The energy equation simplifies to: W_pump = W_potential + W_friction_total
    # This is because P1=P2, v1=0, and v2=0 (based on problem statement)
    W_pump = W_potential + W_friction_total

    print("--- Final Pump Work Calculation ---")
    print("The work of the pump (W_pump) is the sum of the potential energy change and the total friction losses.")
    print("W_pump = g * (z2 - z1) + W_friction_total")
    print(f"W_pump = {W_potential:.4f} J/kg + {W_friction_total:.4f} J/kg")
    print(f"The calculated work of the pump is: {W_pump:.4f} J/kg")

if __name__ == '__main__':
    calculate_pump_work()
    # Final answer needs to be enclosed in <<<>>>
    # After running the code, the final result is ~33.12 J/kg.
    print("\n<<<33.12>>>")