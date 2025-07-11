import math

def solve_b_field():
    """
    Calculates the magnitude of the magnetic field at the point (1, -1, 0)
    due to two infinite current-carrying wires.
    """
    # The problem asks for a symbolic answer in terms of mu_0 (permeability of free space) and I (current).
    # We will perform the calculations for the numerical coefficients.
    
    # Point of interest
    x, y, z = 1, -1, 0
    print(f"Calculating the magnetic field at the point P = ({x}, {y}, {z}).")
    print("-" * 40)

    # --- Wire 1: along x-axis, current in +x direction ---
    print("Analysis for Wire 1 (along x-axis):")
    # Perpendicular distance r1 from the x-axis to P(x, y, z) is sqrt(y^2 + z^2)
    r1 = math.sqrt(y**2 + z**2)
    print(f"The perpendicular distance r1 = sqrt({y}^2 + {z}^2) = {r1}")
    # The magnitude B1 is proportional to 1/r1. The direction is -z.
    b1_z_coeff = -1 / r1
    print(f"The magnetic field vector B1 is proportional to (0, 0, {b1_z_coeff})")
    print("-" * 40)
    
    # --- Wire 2: along y-axis, current in +y direction ---
    print("Analysis for Wire 2 (along y-axis):")
    # Perpendicular distance r2 from the y-axis to P(x, y, z) is sqrt(x^2 + z^2)
    r2 = math.sqrt(x**2 + z**2)
    print(f"The perpendicular distance r2 = sqrt({x}^2 + {z}^2) = {r2}")
    # The magnitude B2 is proportional to 1/r2. The direction is -z.
    b2_z_coeff = -1 / r2
    print(f"The magnetic field vector B2 is proportional to (0, 0, {b2_z_coeff})")
    print("-" * 40)

    # --- Total Magnetic Field ---
    # The total field vector's coefficient is the sum of the individual coefficients,
    # all multiplied by the base unit of (mu_0 * I) / (2 * pi).
    b_total_z_coeff = b1_z_coeff + b2_z_coeff
    
    print("Total Magnetic Field Calculation:")
    print(f"The total magnetic field vector B_total has a z-component coefficient of {b1_z_coeff} + {b2_z_coeff} = {b_total_z_coeff}.")
    print(f"So, B_total_vec = (0, 0, {b_total_z_coeff}) * (mu_0 * I) / (2 * pi)")
    print(f"This simplifies to B_total_vec = (0, 0, {-1}) * (mu_0 * I) / pi\n")
    
    # --- Magnitude of the Total Field ---
    # Magnitude is |b_total_z_coeff| * (mu_0 * I) / (2 * pi)
    magnitude_coeff = abs(b_total_z_coeff)
    
    # Final simplified magnitude is magnitude_coeff/2 * (mu_0 * I) / pi
    numerator_coeff = int(magnitude_coeff / 2)
    denominator_coeff = 1
    
    print("The final simplified equation for the magnitude is:")
    print(f"Magnitude = ({numerator_coeff} * mu_0 * I) / ({denominator_coeff} * pi)")

solve_b_field()