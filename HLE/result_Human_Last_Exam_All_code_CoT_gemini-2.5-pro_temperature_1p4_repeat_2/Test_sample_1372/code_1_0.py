import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """

    # --- Given Data ---
    # Side of the equilateral triangle in meters
    a = 1.2 * 10**10
    # Tangential velocity in km/s
    v_kms = 125
    # The solar mass in kg
    M_solar_kg = 1.99 * 10**30
    # The gravitational constant in kg^-1 m^3 s^-2
    G = 6.67 * 10**-11

    # --- Explanation ---
    # The problem provides instantaneous velocity and orbital period which are inconsistent
    # for a circular orbit. We will proceed by using the instantaneous values (v and a),
    # as the force balance equation holds true at the instant described.

    # The derived formula for the mass of a single component (M) is: M = (v^2 * a) / G
    # To get the answer in solar masses, we divide by the solar mass in kg:
    # M_sol = (v^2 * a) / (G * M_solar_kg)

    # --- Calculations ---
    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # Calculate the numerator of the equation: v^2 * a
    numerator = (v_ms**2) * a

    # Calculate the denominator of the equation: G * M_solar_kg
    denominator = G * M_solar_kg

    # Calculate the final mass in solar masses
    mass_in_solar_masses = numerator / denominator

    # --- Output ---
    print("To find the mass of a single component, we equate the net gravitational force on one star to the centripetal force at the given instant.")
    print("The resulting equation for the mass in solar masses is:")
    print("Mass = (v^2 * a) / (G * M_solar)\n")

    print("Substituting the numerical values into the equation:")
    # The prompt requires printing each number in the final equation.
    # The final equation calculates the mass in solar masses.
    print(f"{mass_in_solar_masses:.1f} solar masses = (({v_kms}*1000 m/s)^2 * {a:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2 * {M_solar_kg:.2e} kg)")

    # The final answer must be returned in the specific format.
    final_answer = round(mass_in_solar_masses, 1)
    
    # This print statement is for display purposes, showing the clean result.
    print(f"\nThe calculated mass of a single component is {final_answer} solar masses.")


solve_star_mass()
print("<<<1.4>>>")