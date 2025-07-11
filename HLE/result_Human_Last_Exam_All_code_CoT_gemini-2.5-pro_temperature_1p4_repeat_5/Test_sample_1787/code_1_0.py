import math

def calculate_brightness_drop():
    """
    Solves for the brightness drop of a brown dwarf during a planetary transit based on a set of observations.

    The solution proceeds in several steps:
    1.  Translate the two observations into mathematical equations.
    2.  Use principles of orbital mechanics to find a relationship between the planets' orbital distances.
    3.  Solve the system of equations for the ratio of the planet's radius to the brown dwarf's radius.
    4.  Calculate the final brightness drop in magnitudes.
    """

    # --- Step 1 & 2: Formulate and analyze the equations ---
    # From Observation 1 (from Space Station/Brown Dwarf):
    # Angular size of P1 = 0.8 * Angular size of P2
    # (2*R1/d1) = 0.8 * (2*R2/d2)  =>  R1/R2 = 0.8 * (d1/d2)

    # From Observation 2 (from Planet 2):
    # Angular size of P1 = 0.2 * Angular size of Brown Dwarf
    # (2*R1/(d2-d1)) = 0.2 * (2*R_bd/d2) => R1/R_bd = 0.2 * (1 - d1/d2)
    #
    # We have two equations but three unknown ratios (R1/R2, R1/R_bd, d1/d2).
    # We need a third relation.

    # --- Step 3: Incorporate Orbital Mechanics ---
    # The problem states P1 is at the pericenter of a parabolic orbit and P2 is in a circular orbit.
    # The pericenter distance (d1) of a parabolic orbit is d1 = h1^2 / (2*G*M).
    # The radius (d2) of a circular orbit is d2 = h2^2 / (G*M).
    # Assuming the specific angular momentums are equal (h1=h2), which is the physical condition
    # that makes the problem solvable, we get:
    # d2 = 2 * d1  =>  d1/d2 = 0.5
    ratio_d1_d2 = 0.5

    # --- Step 4: Solve for the Radius Ratio ---
    # Now we can solve for R2/R_bd.
    # From the equations in Step 1 & 2:
    # R1/R_bd = 0.2 * (1 - ratio_d1_d2) = 0.2 * (1 - 0.5) = 0.1
    # R1/R2 = 0.8 * ratio_d1_d2 = 0.8 * 0.5 = 0.4
    # Therefore, R2/R_bd = (R1/R_bd) / (R1/R2) = 0.1 / 0.4 = 0.25
    ratio_R2_Rbd = 0.25

    # --- Step 5: Calculate Brightness Drop ---
    # The brightness drop depends on the ratio of the disk areas.
    # This is the square of the ratio of the radii.
    ratio_areas = ratio_R2_Rbd ** 2

    # The flux during transit is F_transit = F_initial * (1 - ratio_areas)
    flux_ratio = 1 - ratio_areas

    # The change in magnitude is delta_m = -2.5 * log10(F_transit / F_initial)
    delta_m = -2.5 * math.log10(flux_ratio)

    # The problem asks to output the numbers in the final equation.
    # The final equation is delta_m = -2.5 * log10(1 - (R2/R_bd)^2)
    print(f"The ratio of the squared radii (R2/R_bd)^2, which is the fractional flux drop, is: {ratio_areas}")
    print(f"The brightness drop in bolometric magnitudes is: {delta_m:.3f}")


if __name__ == "__main__":
    calculate_brightness_drop()
    # The final numerical answer.
    # delta_m = -2.5 * math.log10(1 - 0.25**2) = 0.06996...
    # Rounded to three decimal places.
    final_answer = 0.070
    # print(f"<<<{final_answer:.3f}>>>")