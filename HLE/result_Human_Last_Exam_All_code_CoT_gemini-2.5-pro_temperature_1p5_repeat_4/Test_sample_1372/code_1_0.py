import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    side_length_s = 1.2 * 10**10  # m
    velocity_v_km_s = 125          # km/s
    G = 6.67 * 10**-11             # m^3 kg^-1 s^-2
    solar_mass = 1.99 * 10**30     # kg
    
    # --- Step 1: Convert units ---
    velocity_v_m_s = velocity_v_km_s * 1000 # convert km/s to m/s
    
    # --- Step 2: Derive and calculate the mass m ---
    # The gravitational force on one star due to the other two is F_g = sqrt(3) * G * m^2 / s^2.
    # The radius of the orbit around the center of mass is r = s / sqrt(3).
    # The centripetal force is F_c = m * v^2 / r = m * v^2 * sqrt(3) / s.
    # Equating F_g and F_c and simplifying gives: G * m / s = v^2.
    # Solving for m, we get the final formula: m = v^2 * s / G.
    
    mass_kg = (velocity_v_m_s**2 * side_length_s) / G
    
    # --- Step 3: Print the equation with numerical values ---
    print("The formula to calculate the mass (m) of a single component is: m = v^2 * s / G")
    print("\nPlugging in the values:")
    print(f"m = ({velocity_v_m_s:.2e} m/s)^2 * ({side_length_s:.2e} m) / ({G:.2e} m^3 kg^-1 s^-2)")
    
    # --- Step 4: Calculate and print the result in kg and solar masses ---
    print(f"\nCalculated mass in kilograms: {mass_kg:.4e} kg")
    
    mass_in_solar_masses = mass_kg / solar_mass
    
    print("\nConverting to solar masses:")
    print(f"Mass in solar masses = {mass_kg:.4e} kg / {solar_mass:.2e} kg")
    print(f"Result: {mass_in_solar_masses:.4f} solar masses")
    
    # --- Step 5: Print the final answer rounded to one decimal place ---
    print("\n---")
    print(f"The mass of a single component, accurate to one decimal place, is {mass_in_solar_masses:.1f} solar masses.")
    
calculate_star_mass()

print("<<<1.4>>>")