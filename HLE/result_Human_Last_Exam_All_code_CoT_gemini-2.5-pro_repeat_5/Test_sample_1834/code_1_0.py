import math

def calculate_magnetic_field_magnitude():
    """
    This function calculates the magnitude of the magnetic field at the point (1, -1, 0)
    due to two infinite current-carrying wires.
    """
    
    # Point of interest
    px, py, pz = 1, -1, 0
    
    print("The goal is to find the magnitude of the total magnetic field B_total = B1 + B2.")
    print("The magnitude of the magnetic field from a single infinite wire is B = (μ₀ * I) / (2 * π * r).\n")

    # --- Calculation for Wire 1 (on x-axis) ---
    print("Step 1: Analyze the field B1 from Wire 1 (on the x-axis).")
    # The perpendicular distance r1 is the distance from the point to the x-axis.
    r1 = math.sqrt(py**2 + pz**2)
    print(f"The point is P(x,y,z) = ({px}, {py}, {pz}).")
    print(f"The perpendicular distance to the x-axis is r1 = sqrt(y^2 + z^2) = sqrt(({py})^2 + {pz}^2) = {r1}.")
    print("The magnitude of B1 is |B1| = (μ₀ * I) / (2 * π * r1).")
    print("Using the right-hand rule, the direction of B1 at this point is the negative z-direction.\n")

    # --- Calculation for Wire 2 (on y-axis) ---
    print("Step 2: Analyze the field B2 from Wire 2 (on the y-axis).")
    # The perpendicular distance r2 is the distance from the point to the y-axis.
    r2 = math.sqrt(px**2 + pz**2)
    print(f"The point is P(x,y,z) = ({px}, {py}, {pz}).")
    print(f"The perpendicular distance to the y-axis is r2 = sqrt(x^2 + z^2) = sqrt(({px})^2 + {pz}^2) = {r2}.")
    print("The magnitude of B2 is |B2| = (μ₀ * I) / (2 * π * r2).")
    print("Using the right-hand rule, the direction of B2 at this point is also the negative z-direction.\n")

    # --- Summation ---
    print("Step 3: Calculate the magnitude of the total magnetic field.")
    print("Since both B1 and B2 point in the same direction (-z), we can add their magnitudes directly.")
    print("Magnitude |B_total| = |B1| + |B2|")
    print(f"|B_total| = (μ₀ * I) / (2 * π * {r1}) + (μ₀ * I) / (2 * π * {r2})")
    print(f"|B_total| = (μ₀ * I) * (1 / (2 * π * {r1}) + 1 / (2 * π * {r2}))")
    print(f"|B_total| = (μ₀ * I) * (2 / (2 * π))")
    print("|B_total| = (μ₀ * I) / π")

# Execute the function to print the derivation
calculate_magnetic_field_magnitude()