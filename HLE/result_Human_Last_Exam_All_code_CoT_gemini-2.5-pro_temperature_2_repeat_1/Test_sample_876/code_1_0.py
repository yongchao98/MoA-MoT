import math

def calculate_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # --- Given values ---
    r_w = 11.25  # Wire radius in mm
    m = 17.32    # Distance from wire center to cable center in mm
    R = 32.32    # External radius of the cable in mm
    epsilon_r = 4.2  # Relative permittivity of the insulator (dimensionless)

    # --- Physical Constant ---
    epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

    # --- Calculations ---
    # The formula for capacitance per phase is C = (2 * pi * ε₀ * ε_r) / ln((R² - m²) / (R * r_w))
    
    print("This script calculates the capacitance of a three-phase cable with a common screen.")
    print("The formula for capacitance per phase (C) is:")
    print("C = (2 * π * ε₀ * ε_r) / ln((R² - m²) / (R * r_w)) [in F/m]\n")

    print("Step 1: List the given values and constants.")
    print(f"  Wire radius (r_w) = {r_w} mm")
    print(f"  Distance to center (m) = {m} mm")
    print(f"  External radius (R) = {R} mm")
    print(f"  Relative permittivity (ε_r) = {epsilon_r}")
    print(f"  Vacuum permittivity (ε₀) = {epsilon_0} F/m")
    print(f"  Pi (π) ≈ {math.pi}\n")

    # Numerator of the formula
    numerator = 2 * math.pi * epsilon_0 * epsilon_r
    
    # Terms for the denominator's argument
    R_squared = R**2
    m_squared = m**2
    R_times_rw = R * r_w
    log_argument = (R_squared - m_squared) / R_times_rw
    
    # Denominator of the formula
    denominator = math.log(log_argument)
    
    # Capacitance in F/m
    C_per_meter = numerator / denominator
    
    # Capacitance in μF/km
    C_per_km_uF = C_per_meter * 1e9

    print("Step 2: Substitute the values into the equation.\n")
    print("First, calculate the numerator: 2 * π * ε₀ * ε_r")
    print(f"  = 2 * {math.pi:.5f} * {epsilon_0} * {epsilon_r}")
    print(f"  = {numerator}\n")
    
    print("Next, calculate the denominator: ln((R² - m²) / (R * r_w))")
    print(f"  = ln(({R}² - {m}²) / ({R} * {r_w}))")
    print(f"  = ln(({R_squared:.4f} - {m_squared:.4f}) / {R_times_rw:.4f})")
    print(f"  = ln({(R_squared - m_squared):.4f} / {R_times_rw:.4f})")
    print(f"  = ln({log_argument:.5f})")
    print(f"  = {denominator:.5f}\n")

    print("Step 3: Calculate the capacitance in Farads per meter (F/m).")
    print(f"C = {numerator} / {denominator}")
    print(f"C = {C_per_meter} F/m\n")

    print("Step 4: Convert the result to microfarads per kilometer (μF/km).")
    print("  C (μF/km) = C (F/m) * 10^9")
    print(f"  C = {C_per_meter} * 10^9")
    print(f"  C = {C_per_km_uF:.4f} μF/km\n")

    print("--- FINAL RESULT ---")
    print(f"The capacitance of the three-phase cable is {C_per_km_uF:.4f} μF/km.")

if __name__ == "__main__":
    calculate_capacitance()
<<<0.3260>>>