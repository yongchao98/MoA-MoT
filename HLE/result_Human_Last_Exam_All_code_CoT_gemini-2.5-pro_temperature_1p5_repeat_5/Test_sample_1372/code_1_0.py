import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # --- Given Data ---
    # Gravitational constant in kg^-1 m^3 s^-2
    G = 6.67e-11
    # Solar mass in kg
    M_sun = 1.99e30
    # Side of the equilateral triangle in meters
    a = 1.2e10
    # Tangential velocity of each star in m/s (125 km/s)
    v = 125 * 1000

    print("Step 1: Determine the equation for the mass (m) of a single star.")
    print("The centripetal force (Fc = m*v^2/r) on one star is provided by the net gravitational force from the other two stars (Fnet = sqrt(3)*G*m^2/a^2).")
    print("The orbital radius is r = a / sqrt(3).")
    print("By equating Fnet = Fc and simplifying, we get the formula:")
    print("m = (v^2 * a) / G\n")

    print("Step 2: Substitute the given values into the formula to find the mass in kg.")
    # Calculate the mass in kilograms
    mass_kg = (v**2 * a) / G
    
    print(f"m = ({v} m/s)^2 * ({a} m) / ({G} m^3 kg^-1 s^-2)")
    print(f"m = ({v**2}) * ({a}) / ({G})")
    print(f"m = {v**2 * a} / {G}")
    print(f"m = {mass_kg:.4e} kg\n")

    print("Step 3: Convert the mass from kg to solar masses.")
    # Convert mass to solar masses
    mass_solar = mass_kg / M_sun
    print(f"Mass in solar masses = Mass in kg / Solar Mass")
    print(f"Mass in solar masses = {mass_kg:.4e} kg / {M_sun:.4e} kg")
    print(f"Mass in solar masses = {mass_solar:.4f}\n")
    
    # Round the final answer to one decimal place
    final_answer = round(mass_solar, 1)

    print(f"Final Answer: The mass of a single component, rounded to one decimal place, is {final_answer} solar masses.")
    
    # The final required output format
    print(f"<<<{final_answer}>>>")

calculate_star_mass()