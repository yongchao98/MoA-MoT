import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at point (1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Symbolic representations for clarity in output
    mu_0 = "μ₀"
    I = "I"
    pi = "π"

    # Coordinates of the point P
    px, py, pz = 1, -1, 0

    print(f"Calculating the magnetic field at point P({px}, {py}, {pz}).\n")

    # --- Analysis for Wire 1 (on x-axis) ---
    print("--- Step 1: Magnetic Field from Wire 1 (on x-axis) ---")
    # Perpendicular distance r1 from the x-axis to P(1, -1, 0)
    r1 = math.sqrt(py**2 + pz**2)
    print(f"The perpendicular distance r1 from the x-axis to P is sqrt({py}^2 + {pz}^2) = {int(r1)}")
    print(f"The magnitude of the field B1 is given by B1 = ({mu_0} * {I}) / (2 * {pi} * r1).")
    print("By the right-hand rule, with current in the +x direction, the magnetic field at P (a point with a negative y-coordinate) is in the -z direction.")
    print(f"Therefore, the vector B1 = - (({mu_0} * {I}) / (2 * {pi} * {int(r1)})) k\n")

    # --- Analysis for Wire 2 (on y-axis) ---
    print("--- Step 2: Magnetic Field from Wire 2 (on y-axis) ---")
    # Perpendicular distance r2 from the y-axis to P(1, -1, 0)
    r2 = math.sqrt(px**2 + pz**2)
    print(f"The perpendicular distance r2 from the y-axis to P is sqrt({px}^2 + {pz}^2) = {int(r2)}")
    print(f"The magnitude of the field B2 is given by B2 = ({mu_0} * {I}) / (2 * {pi} * r2).")
    print("By the right-hand rule, with current in the +y direction, the magnetic field at P (a point with a positive x-coordinate) is in the -z direction.")
    print(f"Therefore, the vector B2 = - (({mu_0} * {I}) / (2 * {pi} * {int(r2)})) k\n")

    # --- Superposition and Final Calculation ---
    print("--- Step 3: Total Magnetic Field ---")
    print("The total magnetic field B_total is the vector sum of B1 and B2.")
    print("Since both B1 and B2 are in the -z direction, their magnitudes add up.")
    print("Magnitude |B_total| = |B1| + |B2|")
    print("\nThe final equation for the magnitude is:")
    
    # We explicitly show the numbers from the problem being used in the final equation.
    print(f"|B_total| = ({mu_0} * {I}) / (2 * {pi} * {int(r1)}) + ({mu_0} * {I}) / (2 * {pi} * {int(r2)})")
    print(f"|B_total| = 2 * (({mu_0} * {I}) / (2 * {pi}))")
    print(f"|B_total| = ({mu_0} * {I}) / {pi}")

solve_magnetic_field()