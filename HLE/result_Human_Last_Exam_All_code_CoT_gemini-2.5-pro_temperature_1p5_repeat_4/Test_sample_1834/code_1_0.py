import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point due to two infinite wires.
    """
    # Define symbolic constants for the output string
    mu_0_I_div_2pi = "(μ₀ * I) / (2π)"

    # Point of interest P(x, y, z)
    P = (1, -1, 0)
    x, y, z = P

    print("Step-by-step calculation of the magnetic field at P(1, -1, 0):")
    print("="*60)

    # --- Field from Wire 1 (along x-axis, current in +x direction) ---
    print("1. Magnetic Field from Wire 1 (along x-axis):")
    
    # The perpendicular distance from a point (x,y,z) to the x-axis is sqrt(y^2 + z^2)
    r1_val = math.sqrt(y**2 + z**2)
    print(f"   - The point P is (x,y,z) = ({x}, {y}, {z}).")
    print(f"   - The perpendicular distance r1 from P to the x-axis is calculated as sqrt(y² + z²).")
    print(f"   - r1 = sqrt(({y})² + ({z})²) = {r1_val}")
    
    # Magnitude of B1
    B1_mag_str = f"{mu_0_I_div_2pi} / {r1_val}"
    print(f"   - The magnitude of the magnetic field |B1| is given by (μ₀ * I) / (2π * r1).")
    print(f"   - |B1| = {B1_mag_str}")
    
    # Direction of B1 using the right-hand rule
    # Current is in +x. The point is at y=-1 (below the wire). The field points in the -z direction.
    print(f"   - Using the right-hand rule, with current in the +x direction, the magnetic field at a point with a negative y-value points in the -z direction.")
    print(f"   - Therefore, the vector B1 = - |B1| k̂ = -({B1_mag_str}) k̂.")
    print("-"*60)

    # --- Field from Wire 2 (along y-axis, current in +y direction) ---
    print("2. Magnetic Field from Wire 2 (along y-axis):")
    
    # The perpendicular distance from a point (x,y,z) to the y-axis is sqrt(x^2 + z^2)
    r2_val = math.sqrt(x**2 + z**2)
    print(f"   - The perpendicular distance r2 from P to the y-axis is calculated as sqrt(x² + z²).")
    print(f"   - r2 = sqrt(({x})² + ({z})²) = {r2_val}")

    # Magnitude of B2
    B2_mag_str = f"{mu_0_I_div_2pi} / {r2_val}"
    print(f"   - The magnitude of the magnetic field |B2| is given by (μ₀ * I) / (2π * r2).")
    print(f"   - |B2| = {B2_mag_str}")

    # Direction of B2 using the right-hand rule
    # Current is in +y. The point is at x=1 (to the right of the wire). The field points in the +z direction.
    print(f"   - Using the right-hand rule, with current in the +y direction, the magnetic field at a point with a positive x-value points in the +z direction.")
    print(f"   - Therefore, the vector B2 = + |B2| k̂ = +({B2_mag_str}) k̂.")
    print("-"*60)
    
    # --- Total Magnetic Field ---
    print("3. Total Magnetic Field (by Superposition):")
    print("   - The total magnetic field B_total is the vector sum B1 + B2.")
    print(f"   - B_total = [ -({mu_0_I_div_2pi} / {r1_val}) ] k̂  +  [ +({mu_0_I_div_2pi} / {r2_val}) ] k̂")
    
    # The final equation and numbers
    # Since r1 and r2 are both 1, the magnitudes are equal and the vectors are opposite.
    final_magnitude = 0.0
    print(f"   - Since {r1_val} = {r2_val} = 1, the magnitudes are identical, and the directions are opposite. The vectors cancel out.")
    print(f"   - B_total = 0")
    print("="*60)

    # --- Final Answer: Magnitude ---
    print("\nFinal Answer:")
    print(f"The magnitude of the total magnetic field at ({x},{y},{z}) is {final_magnitude}.")

solve_magnetic_field()