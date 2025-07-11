import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
    L = 1.2e10    # Side of the equilateral triangle in meters
    v_kms = 125   # Tangential velocity in km/s
    M_sun = 1.99e30 # Solar mass in kg

    # Convert velocity to m/s
    v_ms = v_kms * 1000

    # The problem provides both velocity (v) and period (T), which are inconsistent
    # for the given geometry. We use the velocity 'v' as it describes the
    # instantaneous state of the system, which directly relates to the forces involved.
    # The governing equation is derived by equating the net gravitational force on a star
    # with the centripetal force required for its orbit.
    # F_gravity = F_centripetal
    # sqrt(3) * G * m^2 / L^2 = m * v^2 / r, where r = L / sqrt(3)
    # This simplifies to the equation: m = v^2 * L / G

    # Calculate the mass of a single component in kg
    m_kg = (v_ms**2 * L) / G

    # Convert the mass to solar masses
    m_solar = m_kg / M_sun

    # Print the explanation and the result
    print("The equation for the mass (m) of a single component is derived from F_gravity = F_centripetal.")
    print("This simplifies to: m = v^2 * L / G")
    print("\nPlugging in the values:")
    print(f"m = ({v_ms:.2e} m/s)^2 * ({L:.2e} m) / ({G:.2e} m^3 kg^-1 s^-2)")
    print(f"m = {m_kg:.4e} kg")
    print(f"\nConverting to solar masses (M_sun = {M_sun:.2e} kg):")
    print(f"Mass in solar masses = {m_kg:.4e} kg / {M_sun:.2e} kg")
    print(f"\nFinal Answer: The mass of a single component is {m_solar:.1f} solar masses.")


calculate_star_mass()
print("\n<<<1.4>>>")