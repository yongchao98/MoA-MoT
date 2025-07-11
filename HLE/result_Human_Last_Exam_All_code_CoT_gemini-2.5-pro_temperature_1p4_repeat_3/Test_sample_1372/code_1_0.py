import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_sun = 1.99e30 # Solar mass in kg
    L = 1.2e10    # Side of the equilateral triangle in m
    v = 125 * 1000 # Tangential velocity in m/s

    # The equation to find the mass 'm' of a single component is derived
    # from equating the net gravitational force with the centripetal force:
    # F_net = F_c
    # (sqrt(3) * G * m^2) / L^2 = m * v^2 / r, where r = L / sqrt(3)
    # This simplifies to: m = (v^2 * L) / G

    # Calculate the mass in kg
    mass_kg = (v**2 * L) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / M_sun

    # Print the calculation steps
    print("Step 1: Define constants and given values.")
    print(f"Gravitational Constant (G): {G} m^3 kg^-1 s^-2")
    print(f"Side of triangle (L): {L:.1e} m")
    print(f"Tangential velocity (v): {v:.3e} m/s")
    print(f"Solar Mass (M_sun): {M_sun:.2e} kg\n")

    print("Step 2: Use the derived formula m = (v^2 * L) / G to find the mass in kg.")
    print(f"m = ({v:.3e}^2 * {L:.1e}) / {G}")
    print(f"m = ({v**2:.3e} * {L:.1e}) / {G}")
    print(f"m = ({(v**2 * L):.3e}) / {G}")
    print(f"m = {mass_kg:.3e} kg\n")

    print("Step 3: Convert the mass from kg to solar masses.")
    print(f"Mass in solar masses = Mass in kg / Solar Mass")
    print(f"Mass in solar masses = {mass_kg:.3e} / {M_sun:.2e}")
    print(f"Mass in solar masses = {mass_solar:.3f}\n")
    
    print("Final Answer:")
    print(f"The mass of a single component is {mass_solar:.1f} solar masses.")

calculate_star_mass()
<<<1.4>>>