import math

def solve_magnetic_field():
    """
    Calculates the magnitude of the magnetic field at a point P(1,-1,0)
    due to two infinite current-carrying wires.
    """
    # --- Problem Parameters ---
    # We use strings for symbolic representation.
    mu_0 = "μ₀"
    I = "I"
    pi = "π"

    # Point P coordinates
    p_x, p_y, p_z = 1, -1, 0

    print("This script calculates the magnitude of the magnetic field at a given point due to two infinite wires.")
    print("---")
    print(f"Point P = ({p_x}, {p_y}, {p_z})")
    print("Wire 1: Placed on the x-axis, current I flows in the +x direction.")
    print("Wire 2: Placed on the y-axis, current I flows in the +y direction.")
    print("---\n")

    # --- Calculation for Wire 1 (on x-axis) ---
    print("Step 1: Calculate the magnetic field from Wire 1 (on x-axis).")
    # Perpendicular distance from P(x,y,z) to the x-axis is sqrt(y^2 + z^2)
    r1 = math.sqrt(p_y**2 + p_z**2)
    print(f"   - The distance r₁ from point P to Wire 1 is √(({p_y})² + ({p_z})²) = {r1}")

    # Direction of B1 using the right-hand rule:
    # Current (I₁) is in +x direction (î). Vector from wire to point P is in -y direction (-ĵ).
    # The direction of B₁ is proportional to î × (-ĵ) = -k̂.
    print("   - Using the right-hand rule, the direction of the magnetic field B₁ is in the negative z-direction (-k̂).")
    print(f"   - The vector B₁ = -({mu_0} * {I}) / (2 * {pi} * {r1}) k̂ = -({mu_0} * {I}) / (2 * {pi}) k̂")
    print("-" * 20)

    # --- Calculation for Wire 2 (on y-axis) ---
    print("Step 2: Calculate the magnetic field from Wire 2 (on y-axis).")
    # Perpendicular distance from P(x,y,z) to the y-axis is sqrt(x^2 + z^2)
    r2 = math.sqrt(p_x**2 + p_z**2)
    print(f"   - The distance r₂ from point P to Wire 2 is √(({p_x})² + ({p_z})²) = {r2}")

    # Direction of B2 using the right-hand rule:
    # Current (I₂) is in +y direction (ĵ). Vector from wire to point P is in +x direction (î).
    # The direction of B₂ is proportional to ĵ × î = -k̂.
    print("   - Using the right-hand rule, the direction of the magnetic field B₂ is also in the negative z-direction (-k̂).")
    print(f"   - The vector B₂ = -({mu_0} * {I}) / (2 * {pi} * {r2}) k̂ = -({mu_0} * {I}) / (2 * {pi}) k̂")
    print("-" * 20)

    # --- Total Magnetic Field ---
    print("Step 3: Calculate the total magnetic field using superposition (B_total = B₁ + B₂).")
    print(f"   - B_total = [ -({mu_0} * {I}) / (2 * {pi}) k̂ ] + [ -({mu_0} * {I}) / (2 * {pi}) k̂ ]")
    # Adding the two vectors
    print(f"   - B_total = -2 * ({mu_0} * {I}) / (2 * {pi}) k̂")
    # Simplifying the expression
    print(f"   - B_total = -({mu_0} * {I}) / ({pi}) k̂")
    print("-" * 20)

    # --- Final Magnitude ---
    print("Step 4: Find the magnitude of the total magnetic field.")
    final_magnitude_str = f"({mu_0} * {I}) / {pi}"
    print(f"   - The magnitude is the absolute value of the vector's coefficient.")
    print(f"   - Final Answer: |B_total| = {final_magnitude_str}")

solve_magnetic_field()