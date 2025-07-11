import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    L = 1.2 * 10**10  # side of the equilateral triangle in meters
    v_km_s = 125  # tangential velocity in km/s
    T_days = 1.9 # period in days (Note: This data is inconsistent with v and L, and will be ignored)
    
    # Constants
    G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
    M_solar = 1.99 * 10**30  # Solar mass in kg

    # Convert velocity to m/s
    v = v_km_s * 1000  # velocity in m/s

    # The net gravitational force on one star is F_net = sqrt(3) * G * m^2 / L^2.
    # The centripetal force is F_c = m * v^2 / r, where r = L / sqrt(3).
    # Equating them (F_net = F_c) and solving for m gives: m = (v^2 * L) / G.
    
    # Calculate the mass of a single component in kg
    m_kg = (v**2 * L) / G

    # Convert the mass to solar masses
    m_solar = m_kg / M_solar

    # Print the explanation and the final equation
    print("The mass of a single star component (m) can be found by equating the net gravitational force to the centripetal force.")
    print("The derived formula for the mass is: m = (v^2 * L) / G")
    print("\nPlugging in the values:")
    print(f"m_kg = ({v:.2e} m/s)^2 * ({L:.2e} m) / ({G:.2e} kg^-1 m^3 s^-2)")
    print(f"m_kg = {m_kg:.4e} kg")
    print("\nConverting to solar masses (M_sun):")
    print(f"m_solar = m_kg / M_sun")
    print(f"m_solar = {m_kg:.4e} kg / {M_solar:.2e} kg")
    print(f"m_solar = {m_solar:.4f}")
    
    # Print the final answer rounded to one decimal place
    print(f"\nThe mass of a single component, accurate to one decimal place, is {m_solar:.1f} solar masses.")
    
    # Final answer in the required format
    print(f"\n<<<{m_solar:.1f}>>>")

calculate_star_mass()