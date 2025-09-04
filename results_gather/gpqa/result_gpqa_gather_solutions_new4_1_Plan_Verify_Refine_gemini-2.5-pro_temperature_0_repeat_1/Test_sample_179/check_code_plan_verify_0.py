import math

def check_energy_calculation():
    """
    Calculates the minimum potential energy of the system and checks it against the provided answer.

    The system consists of:
    - 1 central charge q at point P.
    - 12 charges q on a sphere of radius r, centered at P.
    - q = 2e, r = 2 m.

    The total energy U_total = U_center_sphere + U_sphere_sphere.
    - U_center_sphere is the interaction between the central charge and the 12 on the sphere.
    - U_sphere_sphere is the interaction among the 12 on the sphere. For minimum energy,
      these form an icosahedron (Thomson problem for N=12).

    The formula is U_total = (12 + E_12) * (k * q^2 / r), where E_12 is the
    dimensionless energy constant for the icosahedron.
    """
    # --- Constants (using CODATA 2018 values for high precision) ---
    # Coulomb's constant
    k = 8.9875517923e9  # N·m²/C²
    # Elementary charge
    e = 1.602176634e-19  # C
    # Radius of the sphere
    r = 2.0  # m

    # --- Problem Parameters ---
    # Charge of each particle
    q = 2 * e
    # Number of charges on the sphere
    num_sphere_charges = 12

    # --- Calculation of the Icosahedron's Dimensionless Energy (E_12) ---
    # The energy of the 12 charges on the sphere is U_sphere = E_12 * (k * q^2 / r).
    # E_12 is the sum of inverse distances for all pairs on a unit-radius sphere.
    # An icosahedron has 66 pairs of vertices with 3 distinct distances.
    
    # Golden ratio
    phi = (1 + math.sqrt(5)) / 2

    # For a unit-radius icosahedron (r=1):
    # 1. Distance to 5 adjacent vertices (edge length)
    # A common formula for edge length 's' is s = 4 / sqrt(10 + 2*sqrt(5))
    d1_norm = 4 / math.sqrt(10 + 2 * math.sqrt(5))
    # 2. Distance to 5 next-nearest vertices
    d2_norm = d1_norm * phi
    # 3. Distance to 1 antipodal vertex (diameter)
    d3_norm = 2.0

    # There are 30 pairs at distance d1, 30 at d2, and 6 at d3.
    # E_12 = sum(1/d_ij) for all 66 pairs on a unit sphere.
    E_12 = (30 / d1_norm) + (30 / d2_norm) + (6 / d3_norm)

    # --- Total Energy Calculation ---
    # The base energy unit for one pair at distance r
    energy_unit = (k * q**2) / r
    
    # Total energy is the sum of interactions with the center (12 pairs) and
    # interactions on the sphere (E_12 pairs, normalized to radius r).
    total_energy = (num_sphere_charges + E_12) * energy_unit

    # --- Verification ---
    # The provided answer is D, which is 2.822 x 10^-26 J.
    expected_answer_val = 2.822e-26

    # Check if the calculated value is close to the expected value.
    # We use a relative tolerance to account for floating point arithmetic.
    if math.isclose(total_energy, expected_answer_val, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated energy is {total_energy:.4e} J, "
                f"which does not match the expected answer of {expected_answer_val:.4e} J. "
                f"The dimensionless energy constant for the icosahedron was calculated as E_12 = {E_12:.4f}. "
                f"The total energy formula used was U_total = (12 + E_12) * (k * q^2 / r).")

# Run the check
result = check_energy_calculation()
print(result)