import math

def estimate_weight_change():
    """
    Analyzes the change in weight of a running hourglass based on given parameters.
    The dominant positive effect, the impact force of the landing sand, is calculated.
    """

    # --- Given Parameters (converted to SI units) ---
    d = 0.01      # m (from 1 cm)
    h = 0.02      # m (from 2 cm)
    H = 0.04      # m (from 4 cm)
    rho = 1500.0  # kg/m^3
    t = 60.0      # s (from 1 minute)
    g = 9.8       # m/s^2
    pi = math.pi

    print("--- Analysis of Hourglass Weight Change ---")
    print("This script estimates the change in weight (ΔW) of a running hourglass.")
    print("The primary effect considered is the impact force of the sand landing in the bottom chamber.")
    print("\n--- Parameter Values (SI Units) ---")
    print(f"Chamber diameter, d = {d} m")
    print(f"Sand column height, h = {h} m")
    print(f"Chamber height, H = {H} m")
    print(f"Sand density, rho = {rho} kg/m^3")
    print(f"Total fall time, t = {t} s")
    print(f"Gravity, g = {g} m/s^2")

    # --- Formula Derivation ---
    # Mass flow rate (m_dot) = Total Mass / Total Time = (rho * Volume) / t
    # m_dot = (rho * pi * (d/2)^2 * h) / t = (pi * d**2 * h * rho) / (4 * t)
    #
    # Impact velocity (v_f) is maximized when fall distance is H.
    # v_f = sqrt(2 * g * H)
    #
    # Change in weight (Delta_W) = m_dot * v_f
    m_dot = (pi * d**2 * h * rho) / (4 * t)
    v_f = math.sqrt(2 * g * H)
    Delta_W = m_dot * v_f

    # Total sand weight for comparison
    W_sand = (rho * pi * d**2 * h / 4) * g

    print("\n--- Calculation based on the formula: ΔW = (π*d²*h*ρ / 4t) * sqrt(2gH) ---")
    
    # Printing the equation with numerical values as requested
    print("\nFinal Equation with numerical values:")
    print(f"ΔW = (π * {d}**2 * {h} * {rho} / (4 * {t})) * sqrt(2 * {g} * {H})")

    print("\n--- Results ---")
    print(f"Calculated mass flow rate, ṁ = {m_dot:.3e} kg/s")
    print(f"Calculated max impact velocity, v_f = {v_f:.3f} m/s")
    print(f"Estimated weight change, ΔW = {Delta_W:.3e} N (Heavier)")
    print(f"For comparison, total weight of sand, W_sand = {W_sand:.3e} N")
    print(f"The calculated positive weight change corresponds to option D.")

estimate_weight_change()