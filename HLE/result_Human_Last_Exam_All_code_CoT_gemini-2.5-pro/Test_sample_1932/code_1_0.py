import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on a model
    of perfectly elastic sand collisions.
    """
    # Parameters from the problem (converted to SI units)
    d = 0.01  # diameter in meters (1 cm)
    h = 0.02  # height of sand column in meters (2 cm)
    H = 0.04  # height of glass chamber in meters (4 cm)
    rho = 1500  # density of sand in kg/m^3
    t = 60  # time for sand to fall in seconds (1 minute)
    g = 9.8  # acceleration due to gravity in m/s^2

    print("--- Hourglass Weight Change Analysis ---")
    print("This script analyzes the change in weight of an hourglass while it's running.")
    print("We consider the 'largest possible effect', which we model as perfectly elastic collisions of sand grains.\n")

    # Step 1: Calculate the mass flow rate (m_dot)
    sand_volume = (math.pi * d**2 / 4) * h
    sand_mass = rho * sand_volume
    m_dot = sand_mass / t
    
    print("--- Parameters ---")
    print(f"Diameter of sand column (d): {d} m")
    print(f"Height of sand column (h): {h} m")
    print(f"Height of chamber (H): {H} m")
    print(f"Density of sand (rho): {rho} kg/m^3")
    print(f"Total run time (t): {t} s")
    print(f"Gravity (g): {g} m/s^2\n")

    print(f"Total mass of sand: {sand_mass:.4f} kg")
    print(f"Mass flow rate (m_dot): {m_dot:.2e} kg/s\n")
    
    # Step 2: Define the characteristic fall distance for the 'largest effect'
    L = H
    print(f"Characteristic fall distance (L = H): {L} m\n")

    # Step 3: Calculate the impact velocity
    v_impact = math.sqrt(2 * g * L)
    print(f"Impact velocity of sand (v_impact = sqrt(2*g*L)): {v_impact:.3f} m/s\n")
    
    # Step 4: Analyze forces for different models
    # Weight reduction from sand in flight
    W_in_flight = m_dot * v_impact
    # Impact force for inelastic collisions
    F_inelastic = m_dot * v_impact
    # Impact force for elastic collisions
    F_elastic = 2 * m_dot * v_impact
    
    print("--- Force Analysis ---")
    print(f"Weight reduction due to sand in flight (W_in_flight): {W_in_flight:.2e} N")
    print(f"Impact force (inelastic model, F_inelastic): {F_inelastic:.2e} N")
    print(f"Impact force (elastic model, F_elastic): {F_elastic:.2e} N\n")
    
    # Step 5: Calculate the final weight change
    delta_W_inelastic = F_inelastic - W_in_flight
    delta_W_elastic = F_elastic - W_in_flight
    
    print("--- Weight Change (ΔW = F_impact - W_in_flight) ---")
    print(f"ΔW for inelastic model (ideal case): {delta_W_inelastic:.2e} N (Effectively zero)")
    print(f"ΔW for elastic model ('largest effect'): {delta_W_elastic:.2e} N (The hourglass is heavier)\n")
    
    # Step 6: Show the final formula and its components
    print("--- Final Formula Derivation ---")
    print("The estimated weight change ΔW is given by the formula for the elastic model:")
    print("ΔW = m_dot * sqrt(2 * g * L)")
    print("Substituting m_dot = (pi * d^2 * h * rho) / (4 * t) and L = H, we get:")
    print("ΔW = (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * H)")
    
    # Print the formula with numerical values
    # The term (pi * d^2 * h * rho) / (4 * t) is m_dot
    # The term sqrt(2 * g * H) is v_impact
    print("\n--- Final Equation with Values ---")
    print(f"ΔW = (π * ({d})^2 * {h} * {rho} / (4 * {t})) * √(2 * {g} * {H})")
    print(f"ΔW = ({m_dot:.2e}) * √({2 * g * H:.3f})")
    print(f"ΔW = {m_dot:.2e} * {v_impact:.3f}")
    print(f"ΔW = {delta_W_elastic:.2e} N")

if __name__ == '__main__':
    calculate_weight_change()