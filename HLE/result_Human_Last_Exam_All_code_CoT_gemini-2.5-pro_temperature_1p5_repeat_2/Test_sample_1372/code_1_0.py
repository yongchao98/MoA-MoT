import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple-star system.
    """
    # Given data
    # Side of the equilateral triangle in meters
    a = 1.2 * 10**10
    # Tangential velocity in km/s
    v_kms = 125
    
    # Constants
    # Gravitational constant in m^3 kg^-1 s^-2
    G = 6.67 * 10**-11
    # Solar mass in kg
    M_solar = 1.99 * 10**30

    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # Explain the derivation
    print("The mass 'm' of a single component star can be found by equating the net gravitational force to the centripetal force.")
    print("The simplified equation for the mass is: m = (v^2 * a) / G")
    print("\nNote: The provided tangential velocity (125 km/s) and period (1.9 days) are physically inconsistent. The calculation proceeds using the given velocity.")

    # Calculate the mass of a single component in kg using the derived formula
    m_kg = (v_ms**2 * a) / G
    
    # Convert the mass to solar masses
    m_solar = m_kg / M_solar
    
    # Round the final answer to one decimal place
    final_answer = round(m_solar, 1)

    # Print the equation with the numbers plugged in
    print("\nCalculating the mass in kg:")
    print(f"m = ({v_ms:.3e} m/s)^2 * ({a:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2)")
    print(f"m = {m_kg:.4e} kg")

    # Print the conversion to solar masses
    print("\nConverting to solar masses:")
    print(f"Mass in solar masses = {m_kg:.4e} kg / {M_solar:.2e} kg")
    print(f"Mass = {final_answer} solar masses")

    # Output the final answer in the required format
    print(f"\n<<<{final_answer}>>>")

calculate_star_mass()