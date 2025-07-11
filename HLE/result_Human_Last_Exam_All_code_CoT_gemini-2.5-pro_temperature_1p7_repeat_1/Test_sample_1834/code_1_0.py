import numpy as np

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P(x,y,z)
    due to two infinite current-carrying wires.
    """
    # This problem requires us to find the magnetic field at a point P(x,y,z)
    # from two infinite wires. The total magnetic field is the vector sum of the
    # fields from each wire.

    # Let B_unit = (mu_0 * I) / (2 * pi), where mu_0 is the permeability of free space
    # and I is the current. The magnitude of the magnetic field from an infinite
    # wire is B = B_unit / r, where r is the perpendicular distance from the wire.

    # Point of interest
    P = np.array([1, -1, 0])

    # --- Contribution from Wire 1 (along x-axis, current in +x) ---

    # The perpendicular distance (r1) from P to the x-axis is its
    # distance in the y-z plane, which is sqrt(y^2 + z^2).
    r1 = np.sqrt(P[1]**2 + P[2]**2)

    # The magnitude of the magnetic field B1 is B_unit / r1.
    # Its value as a multiple of B_unit is:
    B1_coeff = 1 / r1

    # The direction of B1 is found using the right-hand rule.
    # With current in +x, at a point with y=-1, z=0, the field is in the +z direction.
    B1_direction = np.array([0, 0, 1])

    # The magnetic field vector B1, in units of B_unit.
    B1_vec_coeff = B1_coeff * B1_direction

    # --- Contribution from Wire 2 (along y-axis, current in +y) ---

    # The perpendicular distance (r2) from P to the y-axis is its
    # distance in the x-z plane, which is sqrt(x^2 + z^2).
    r2 = np.sqrt(P[0]**2 + P[2]**2)

    # The magnitude of the magnetic field B2 is B_unit / r2.
    # Its value as a multiple of B_unit is:
    B2_coeff = 1 / r2

    # The direction of B2 is found using the right-hand rule.
    # With current in +y, at a point with x=1, z=0, the field is in the -z direction.
    B2_direction = np.array([0, 0, -1])

    # The magnetic field vector B2, in units of B_unit.
    B2_vec_coeff = B2_coeff * B2_direction

    # --- Total Magnetic Field ---

    # The total magnetic field vector is the sum of the individual vectors.
    B_total_vec_coeff = B1_vec_coeff + B2_vec_coeff

    # The magnitude of the total magnetic field is the norm of the total vector.
    B_total_magnitude_coeff = np.linalg.norm(B_total_vec_coeff)

    # --- Output the results ---

    print("This script calculates the magnitude of the magnetic field at point P(1, -1, 0).")
    print("Let B_unit = (μ₀ * I) / (2π). All field values are expressed as multiples of B_unit.")
    print("-" * 50)

    print("Field from Wire 1 (on x-axis):")
    print(f"Perpendicular distance r1 from P to the x-axis is {r1:.2f}")
    print(f"Magnitude B1 = (1 / {r1:.2f}) * B_unit = {B1_coeff:.2f} * B_unit")
    print(f"Direction is in +z. Vector B1 = {B1_vec_coeff} * B_unit")
    print("-" * 50)

    print("Field from Wire 2 (on y-axis):")
    print(f"Perpendicular distance r2 from P to the y-axis is {r2:.2f}")
    print(f"Magnitude B2 = (1 / {r2:.2f}) * B_unit = {B2_coeff:.2f} * B_unit")
    print(f"Direction is in -z. Vector B2 = {B2_vec_coeff} * B_unit")
    print("-" * 50)

    print("Total Magnetic Field:")
    # The final equation is the vector sum, where only the z-component is non-zero initially.
    print("The total field is the vector sum: B_total = B1 + B2.")
    print("The x and y components are zero for both vectors.")
    print("The equation for the z-component is: B_total_z = B1_z + B2_z")
    print(f"Expressed as multiples of B_unit: {B1_vec_coeff[2]:.2f} + ({B2_vec_coeff[2]:.2f}) = {B_total_vec_coeff[2]:.2f}")
    print("-" * 50)
    
    print(f"The magnitude of the total magnetic field is {B_total_magnitude_coeff:.1f} times (μ₀ * I / (2π)).")
    print(f"Since the two field vectors cancel each other out, the resulting magnitude is exactly 0.")

solve_magnetic_field()