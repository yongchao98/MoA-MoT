import math

def calculate_incident_power():
    """
    Calculates the total power incident on the photovoltaic cell at point B.
    """

    # Data given in the problem
    P_source = 1e9  # Luminosity of the source, P = 1 GW = 10^9 W
    S_cell = 10.0      # Area of the photovoltaic cell, S = 10 m^2
    A_mirror = 1e6     # Area of each satellite's mirror, 1 km^2 = 10^6 m^2
    M_moon = 7.35e22   # Mass of the Moon in kg
    R_moon = 1738e3    # Radius of the Moon in m (1738 km)
    T_orbit = 12 * 3600  # Orbital period in s (12 hours)
    G = 6.67e-11       # Gravitational constant in m^3 kg^-1 s^-2
    pi = math.pi

    print("Step 1: Calculate the semi-major axis 'a' of the satellite orbit.")
    print("Using Kepler's Third Law: T^2 / a^3 = 4 * pi^2 / (G * M)")
    a_cubed = (G * M_moon * T_orbit**2) / (4 * pi**2)
    a = a_cubed**(1/3)
    print(f"The calculated semi-major axis a = {a:,.2f} m.\n")

    print("Step 2: Analyze the light path and derive the power formula.")
    print("The geometry (A-X zenith, Y-B zenith) places X and Y at opposite apsides (apoapsis and periapsis).")
    print("Light from point source A is reflected by mirror X, creating virtual source A'.")
    print("Light from A' is reflected by mirror Y, creating virtual source A''.")
    print("The power P' on cell B is found by tracing the intensity of the diverging beam.")
    print("Following this model, the power P' is given by the formula:")
    print("P' = (P * S) / (8 * pi * (2*a - R)^2)")
    print("Notably, this result is independent of the mirror area and the orbit's eccentricity.\n")
    
    print("Step 3: Calculate the final power P'.")
    # This is the derived formula P' = P*S / (8*pi*(2a-R)^2)
    
    # Calculate each term of the equation.
    numerator = P_source * S_cell
    # Using 'a' and 'R_moon' calculated/given above.
    term_in_parentheses = 2 * a - R_moon
    denominator = 8 * pi * term_in_parentheses**2
    
    P_prime = numerator / denominator

    # Convert the result to microwatts (1 W = 1,000,000 uW)
    P_prime_microwatts = P_prime * 1e6

    print(f"Substituting the values into the formula:")
    print(f"P' = ({P_source:.1e} W * {S_cell:.1f} m^2) / (8 * pi * (2 * {a:,.0f} m - {R_moon:,.0f} m)^2)")
    print(f"P' = {numerator:.1e} / (8 * pi * ({term_in_parentheses:,.0f})^2)")
    print(f"P' = {P_prime:.4e} W")
    print(f"P' = {P_prime_microwatts:.1f} microwatts.\n")
    
    # Final answer rounded to one decimal place.
    final_answer = round(P_prime_microwatts, 1)
    print(f"The total power P' incident on the cell is {final_answer} microwatts.")

    return final_answer

# Execute the calculation and print the final result
final_power = calculate_incident_power()
# The final answer in the requested format
# print(f"<<<{final_power}>>>")