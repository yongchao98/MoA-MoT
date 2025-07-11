import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # Given data and constants
    G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
    solar_mass_kg = 1.99 * 10**30  # Mass of the sun in kg
    a = 1.2 * 10**10  # Side of the equilateral triangle in meters
    v_kms = 125       # Tangential velocity in km/s

    # Convert velocity to m/s for SI unit consistency
    v_ms = v_kms * 1000

    # The formula to find the mass of a single component is m = v^2 * a / G.
    
    # Calculate the mass in kilograms
    mass_kg = (v_ms**2 * a) / G
    
    # Convert the mass to solar masses
    mass_solar = mass_kg / solar_mass_kg
    
    # --- Output ---
    print("The equation to calculate the mass (m) of a single star is:")
    print("m = v^2 * a / G\n")
    
    print("Plugging in the numbers:")
    print(f"m = ({v_ms:.0f} m/s)^2 * ({a:.1e} m) / ({G:.2e} kg^-1 m^3 s^-2)")
    
    # Showing the calculation step
    v_squared = v_ms**2
    v_squared_times_a = v_squared * a
    print(f"m = ({v_squared:.4e}) * ({a:.1e}) / ({G:.2e})")
    print(f"m = ({v_squared_times_a:.4e}) / ({G:.2e})")
    print(f"m = {mass_kg:.4e} kg\n")

    print("Converting the mass to solar masses:")
    print(f"Mass in solar masses = {mass_kg:.4e} kg / {solar_mass_kg:.2e} kg")
    print(f"Mass in solar masses = {mass_solar:.3f}\n")
    
    print("The final answer, accurate to one decimal place:")
    print(f"The mass of a single component is {mass_solar:.1f} solar masses.")

# Run the calculation
calculate_star_mass()