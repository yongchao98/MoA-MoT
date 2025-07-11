import math

def calculate_magnetic_field():
    """
    Calculates and explains the magnitude of the magnetic field at point P(1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # Point of interest
    x, y, z = 1, -1, 0

    # --- Wire 1 (on x-axis, current in +x) ---
    # The perpendicular distance r1 from the x-axis to P(x,y,z) is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)

    # --- Wire 2 (on y-axis, current in +y) ---
    # The perpendicular distance r2 from the y-axis to P(x,y,z) is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)

    print("Step 1: Calculate the magnetic field from each wire.")
    print("The formula for the magnitude of the magnetic field from an infinite wire is B = (μ₀ * I) / (2 * π * r).")
    print(f"For Wire 1 (on x-axis), the perpendicular distance r1 to the point (1, -1, 0) is {r1}.")
    print(f"For Wire 2 (on y-axis), the perpendicular distance r2 to the point (1, -1, 0) is {r2}.")
    print("The magnetic fields from both wires point in the same direction (-z). Therefore, their magnitudes add up.")
    print("\nStep 2: Sum the magnitudes to find the total magnitude.")
    print("Total Magnitude B_total = B_wire1 + B_wire2")
    print("B_total = (μ₀ * I) / (2 * π * r1) + (μ₀ * I) / (2 * π * r2)")
    
    # Printing the equation with the calculated numbers
    print("\nStep 3: Substitute the distance values into the equation.")
    print(f"B_total = (μ₀ * I) / (2 * π * {r1}) + (μ₀ * I) / (2 * π * {r2})")
    
    # Simplify the expression
    print("\nStep 4: Simplify the expression.")
    print(f"Since the denominators are the same, we can combine the terms:")
    print(f"B_total = 2 * (μ₀ * I) / (2 * π * {r1})")
    print("The '2's in the numerator and denominator cancel out.")

    print("\nFinal Answer:")
    print("The final expression for the magnitude of the magnetic field is:")
    print(f"B_total = (μ₀ * I) / (π * {r1})")
    print("Or more simply:")
    print("B_total = (μ₀ * I) / π")

calculate_magnetic_field()