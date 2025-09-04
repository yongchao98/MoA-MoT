import numpy as np
import scipy.constants as const

def check_energy_calculation():
    """
    Checks the correctness of the calculated minimum energy for the given system of charges.
    The code recalculates the value from first principles and compares it to the provided answer.
    """
    # --- Define Physical Constants and Problem Parameters ---
    # Use high-precision values from scipy.constants for accuracy.
    k = 1 / (4 * np.pi * const.epsilon_0)  # Coulomb's constant in N m^2 / C^2
    e = const.e                           # Elementary charge in C
    r = 2.0                               # Radius of the sphere in meters
    q = 2 * e                             # Charge of each particle in C

    # --- Calculate the Dimensionless Energy Constant E_12 for the Icosahedron ---
    # The minimum energy configuration for 12 charges on a sphere (Thomson problem)
    # is an icosahedron. The potential energy of this shell is U_sphere = E_12 * (k * q^2 / r).
    # E_12 can be calculated from the geometry of an icosahedron inscribed in a unit sphere.
    
    # Golden ratio
    phi = (1 + np.sqrt(5)) / 2
    
    # The edge length 's_norm' of an icosahedron inscribed in a unit radius sphere.
    s_norm = 4 / np.sqrt(10 + 2 * np.sqrt(5))
    
    # E_12 is the sum of inverse distances for all 66 pairs of vertices on the unit sphere.
    # These pairs fall into three groups by distance:
    # - 30 pairs at distance s_norm (edges)
    # - 30 pairs at distance s_norm * phi (next-nearest neighbors)
    # - 6 pairs at distance 2.0 (antipodal vertices, i.e., the diameter)
    E_12 = (30 / s_norm) + (30 / (s_norm * phi)) + (6 / 2.0)

    # --- Calculate the Total Minimum Energy ---
    # The total energy is the sum of the energy of the sphere charges interacting with each other
    # and the energy of the sphere charges interacting with the central charge.
    # U_total = U_sphere-sphere + U_center-sphere
    # U_total = (E_12 * k * q**2 / r) + (12 * k * q**2 / r)
    # U_total = (12 + E_12) * (k * q**2 / r)
    
    calculated_energy = (12 + E_12) * (k * q**2) / r

    # --- Compare with the Provided Answer ---
    # The LLM's final answer is <<<A>>>, which corresponds to the value 2.822 x 10^-26 J.
    answer_value = 2.822e-26

    # The problem asks for the answer correct to three decimals. This implies rounding
    # the final value when presented in scientific notation.
    # We check if our calculated value, when formatted to three decimal places, matches the answer.
    
    # Format the calculated value to the same precision as the answer options (X.XXX e-YY).
    formatted_calculated_value = float(f"{calculated_energy:.3e}")

    if formatted_calculated_value == answer_value:
        return "Correct"
    else:
        reason = (
            f"The calculation is incorrect. "
            f"The calculated energy is {calculated_energy:.4e} J. "
            f"When formatted to three decimal places (as per the options), this is {formatted_calculated_value:.3e} J. "
            f"The provided answer is option A, which has a value of {answer_value:.3e} J. "
            f"The calculated value does not match the answer's value."
        )
        return reason

# Execute the check and print the result.
print(check_energy_calculation())