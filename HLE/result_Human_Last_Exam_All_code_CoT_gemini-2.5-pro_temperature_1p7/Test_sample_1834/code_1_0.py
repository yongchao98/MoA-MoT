import math

def solve_magnetic_field():
    """
    This function calculates the magnitude of the magnetic field from two perpendicular
    infinite wires and explains the steps.
    """
    
    # The problem is symbolic, so we will use strings to represent the formula components.
    # mu_0 is the permeability of free space.
    # I is the current.
    
    point_p = (1, -1, 0)
    
    print("This script calculates the magnitude of the magnetic field at point P(1, -1, 0).")
    print("---")

    # Step 1: Analyze Wire 1 (on x-axis)
    print("Step 1: Calculate the magnetic field from Wire 1 (on the x-axis).")
    r1 = math.sqrt(point_p[1]**2 + point_p[2]**2)
    print(f"The perpendicular distance (r1) from P{point_p} to the x-axis is sqrt((-1)^2 + 0^2) = {r1}.")
    print("The magnitude of the field from an infinite wire is B = (mu_0 * I) / (2 * pi * r).")
    print(f"So, the magnitude B1 = (mu_0 * I) / (2 * pi * {int(r1)}) = (mu_0 * I) / (2 * pi).")
    print("Using the right-hand rule (current in +x direction), the magnetic field B1 at P points in the +z direction.")
    print("Therefore, the vector B1 = (0, 0, mu_0*I / (2*pi)).")
    print("---")

    # Step 2: Analyze Wire 2 (on y-axis)
    print("Step 2: Calculate the magnetic field from Wire 2 (on the y-axis).")
    r2 = math.sqrt(point_p[0]**2 + point_p[2]**2)
    print(f"The perpendicular distance (r2) from P{point_p} to the y-axis is sqrt(1^2 + 0^2) = {r2}.")
    print(f"The magnitude B2 = (mu_0 * I) / (2 * pi * {int(r2)}) = (mu_0 * I) / (2 * pi).")
    print("Using the right-hand rule (current in +y direction), the magnetic field B2 at P also points in the +z direction.")
    print("Therefore, the vector B2 = (0, 0, mu_0*I / (2*pi)).")
    print("---")
    
    # Step 3: Apply Superposition Principle
    print("Step 3: Calculate the total magnetic field by vector addition (superposition).")
    print("B_total = B1 + B2")
    print("B_total = (0, 0, mu_0*I / (2*pi)) + (0, 0, mu_0*I / (2*pi))")
    print("B_total = (0, 0, 2 * mu_0*I / (2*pi))")
    print("B_total = (0, 0, mu_0*I / pi)")
    print("---")
    
    # Step 4: Final Magnitude
    print("Step 4: Find the magnitude of the total magnetic field.")
    print("The magnitude |B_total| is the length of the B_total vector.")
    print("Final equation for the magnitude: |B_total| = mu_0 * I / pi")
    print("\nAs requested, here are the numbers in the final equation:")
    print("The coefficient in the numerator is 1.")
    print("The coefficient of pi in the denominator is 1.")

solve_magnetic_field()