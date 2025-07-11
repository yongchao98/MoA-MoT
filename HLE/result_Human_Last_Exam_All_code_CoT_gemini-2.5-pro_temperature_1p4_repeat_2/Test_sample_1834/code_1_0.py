import math

def calculate_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # The problem provides symbolic current I. The result will be in terms of I.
    # The point of interest is P(x, y, z) = (1, -1, 0).
    point_x, point_y = 1.0, -1.0

    # Physical constant mu_0 (permeability of free space in T*m/A)
    mu_0 = 4 * math.pi * 1e-7

    # --- Field from Wire 1 (on x-axis) ---
    # The perpendicular distance r1 from the x-axis to P(1, -1, 0) is |point_y|.
    r1 = abs(point_y)
    # The direction of the current is +x. By the right-hand rule, the magnetic
    # field B1 at P is in the -z direction.
    # The magnitude B1 is (mu_0 * I) / (2 * pi * r1).

    # --- Field from Wire 2 (on y-axis) ---
    # The perpendicular distance r2 from the y-axis to P(1, -1, 0) is |point_x|.
    r2 = abs(point_x)
    # The direction of the current is +y. By the right-hand rule, the magnetic
    # field B2 at P is also in the -z direction.
    # The magnitude B2 is (mu_0 * I) / (2 * pi * r2).

    # --- Total Magnetic Field ---
    # Since both B1 and B2 are in the same direction (-z), the magnitude of the
    # total field is the sum of the individual magnitudes.
    # B_total = B1 + B2 = (mu_0 * I) / (2 * pi * r1) + (mu_0 * I) / (2 * pi * r2)

    # We can calculate the numerical coefficient for I.
    coefficient = (mu_0 / (2 * math.pi * r1)) + (mu_0 / (2 * math.pi * r2))

    print("The magnitude of the magnetic field B is given by the sum of the contributions from each wire.")
    print("Final Equation: B_magnitude = (mu_0 * I) / (2 * pi * r1) + (mu_0 * I) / (2 * pi * r2)")
    print(f"Substituting the distances r1 = {r1} m and r2 = {r2} m:")
    print(f"B_magnitude = (mu_0 * I) / (2 * pi * {r1}) + (mu_0 * I) / (2 * pi * {r2})")
    print("This simplifies to:")
    print("B_magnitude = mu_0 * I / pi")
    print("\nThe numerical value for the coefficient (mu_0 / pi) is:")
    print(f"B_magnitude = ({coefficient:.1e}) * I Tesla")

calculate_magnetic_field()
