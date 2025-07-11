import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # --- Given Data ---
    # Side of the equilateral triangle (m)
    L = 1.2 * 10**10
    # Tangential velocity (km/s)
    v_kms = 125
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67 * 10**-11
    # Solar mass (kg)
    M_sun = 1.99 * 10**30

    # --- Calculation ---
    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # The force balance equation simplifies to: m = v^2 * L / G
    # We will now calculate the mass 'm' using this equation.

    # Calculate the mass of a single component in kg
    m_kg = (v_ms**2 * L) / G

    # Convert the mass to solar masses
    m_solar = m_kg / M_sun
    
    # --- Output ---
    print("The derived equation for the mass 'm' of a single component is:")
    print("m = v^2 * L / G\n")
    
    print("Substituting the numerical values into the equation:")
    print(f"m = ({v_ms:.2e} m/s)^2 * ({L:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2)\n")

    print(f"The calculated mass of a single component is {m_kg:.3e} kg.")
    print(f"This is equivalent to {m_solar:.3f} solar masses.\n")
    
    print(f"The final answer, rounded to one decimal place, is {m_solar:.1f} solar masses.")

solve_star_mass()