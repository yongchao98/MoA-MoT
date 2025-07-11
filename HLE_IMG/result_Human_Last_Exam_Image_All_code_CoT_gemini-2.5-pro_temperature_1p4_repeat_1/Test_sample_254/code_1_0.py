import math

# Define the constants given in the problem in SI units.
P_source = 1e9  # Power of the source: 1 GW = 10^9 W
M_moon = 7.35e22  # Mass of the Moon: 7.35 * 10^22 kg
R_moon = 1738e3  # Radius of the Moon: 1738 km = 1.738 * 10^6 m
G = 6.67e-11  # Gravitational constant: 6.67 * 10^-11 m^3 kg^-1 s^-2
T_orbit = 12 * 3600  # Orbital period: 12 hours = 43200 s
S_cell = 10  # Area of the photovoltaic cell: 10 m^2

def solve_power_on_cell():
    """
    This function calculates the total power incident on the photovoltaic cell at point B.
    """

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # T^2 = (4 * pi^2 * a^3) / (G * M)
    # Rearranging for 'a': a = ( (G * M * T^2) / (4 * pi^2) )^(1/3)
    # For the problem to be solvable without the eccentricity, we assume a circular orbit,
    # where the semi-major axis 'a' is equal to the constant orbital radius.
    a = ((G * M_moon * T_orbit**2) / (4 * math.pi**2))**(1/3)

    # Step 2: Derive and apply the formula for the final power P'.
    # Based on a beam propagation model for the light path A -> X -> Y -> B,
    # where the light spreads from virtual sources created by reflections,
    # a simplified formula for the final power P' on the cell at B is derived.
    # The derivation involves significant cancellations of geometric terms.
    # The final formula is: P' = (P_source * S_cell) / (8 * pi * (2*a - R_moon)^2)
    dist_term = 2 * a - R_moon
    P_final_watts = (P_source * S_cell) / (8 * math.pi * dist_term**2)

    # Step 3: Convert the final power to microwatts for the answer.
    P_final_microwatts = P_final_watts * 1e6

    # Step 4: Print the calculation steps and the final answer.
    print("Step-by-step calculation:")
    print(f"1. Gravitational Constant (G): {G:.2e} m^3 kg^-1 s^-2")
    print(f"2. Mass of the Moon (M): {M_moon:.2e} kg")
    print(f"3. Radius of the Moon (R): {R_moon:.0f} m")
    print(f"4. Orbital Period (T): {T_orbit:.0f} s")
    print(f"5. From Kepler's Law, the calculated semi-major axis (a) is: {a:.1f} m")

    print("\nFinal power calculation formula derived from the beam propagation model:")
    print("P' = (P_source * S_cell) / (8 * pi * (2*a - R)^2)")

    print("\nSubstituting the numerical values into the equation:")
    equation_str = (
        f"P' = ({P_source:.1e} W * {S_cell} m^2) / "
        f"(8 * {math.pi:.5f} * (2 * {a:.1f} m - {R_moon:.0f} m)^2)"
    )
    print(equation_str)

    # Display the evaluation of the formula parts
    evaluated_str = (
        f"P' = {P_source * S_cell:.1e} / (8 * {math.pi:.5f} * ({dist_term:.1f})^2)"
    )
    print(evaluated_str)

    evaluated_str_2 = (
        f"P' = {P_source * S_cell:.1e} / ({8 * math.pi * dist_term**2:.3e})"
    )
    print(evaluated_str_2)
    print(f"P' = {P_final_watts:.3e} W")

    print("\nFinal Answer:")
    print(f"The total power incident on the cell is {P_final_microwatts:.1f} microwatts.")
    return P_final_microwatts

# Execute the function to get the solution.
final_answer = solve_power_on_cell()
print(f"<<<{final_answer:.1f}>>>")