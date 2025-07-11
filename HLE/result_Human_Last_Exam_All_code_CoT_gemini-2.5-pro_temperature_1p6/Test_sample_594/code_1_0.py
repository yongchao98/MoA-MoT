import math

def calculate_stable_pore_radius():
    """
    Calculates the stable radius of a pore with trapped gas during sintering.

    During sintering, the driving force for pore shrinkage (densification) is
    the reduction of surface energy. This creates a sintering pressure (Ps)
    inversely proportional to the pore radius (r): Ps = 2*γ/r, where γ is
    the surface energy.

    If a coarsening gas is trapped in the pore, it exerts a counteracting
    pressure (Pgas) according to the ideal gas law: Pgas = nRT/V, where V is
    the pore volume (4/3)*π*r³.

    Densification halts when Ps = Pgas. This script calculates the equilibrium
    pore radius where this condition is met, demonstrating the physical mechanism
    that leads to effects like lower densities and voids (Options A, B, C, E, F).

    The effect described in option (D) is unlikely because these stable pores in the
    interior would pin grain boundaries, leading to SMALLER grains, not larger ones,
    compared to the surface where gas can escape.
    """
    # --- Input Parameters ---
    # Moles of trapped gas (a hypothetical small amount in one pore)
    n = 1.0e-15  # moles
    # Ideal gas constant
    R = 8.314  # J/(mol·K)
    # Sintering temperature for a ceramic oxide like Alumina (1600 C)
    T = 1600 + 273.15  # Kelvin
    # Surface energy of Alumina
    gamma = 1.0  # J/m^2

    # --- Calculation ---
    # The equilibrium condition is 2*γ/r = nRT / ((4/3)πr³)
    # Solving for r: r = sqrt(3*n*R*T / (8*π*γ))
    numerator = 3 * n * R * T
    denominator = 8 * math.pi * gamma
    stable_radius_squared = numerator / denominator
    stable_radius_m = math.sqrt(stable_radius_squared)
    
    # Convert meters to micrometers for easier interpretation
    stable_radius_microns = stable_radius_m * 1e6

    # --- Output ---
    print("This script calculates the stable pore radius due to trapped gas.")
    print("-" * 60)
    print("EQUILIBRIUM EQUATION:")
    print("Sintering Pressure (2*γ/r) = Gas Pressure (nRT/V)")
    print("\nCALCULATION FOR STABLE RADIUS (r):")
    print("r = sqrt( (3 * n * R * T) / (8 * π * γ) )")
    
    # Output the equation with the final numbers plugged in
    print("\nPLUGGING IN THE VALUES:")
    # Formatting numbers for clarity
    n_str = f"{n:.1e}"
    R_str = f"{R:.3f}"
    T_str = f"{T:.2f}"
    gamma_str = f"{gamma:.1f}"
    pi_str = f"{math.pi:.4f}"
    
    print(f"r = sqrt( (3 * {n_str} * {R_str} * {T_str}) / (8 * {pi_str} * {gamma_str}) )")

    # Final result
    print("-" * 60)
    print(f"The calculated stable pore radius is {stable_radius_microns:.2f} micrometers.")
    print("This demonstrates that trapped gas can prevent pores from closing, leading to most of the listed effects.")

calculate_stable_pore_radius()