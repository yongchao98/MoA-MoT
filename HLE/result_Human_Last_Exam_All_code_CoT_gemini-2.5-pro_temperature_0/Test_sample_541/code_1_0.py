import math

def solve_isoperimetric_problem():
    """
    This function demonstrates the isoperimetric problem by comparing the areas
    of different shapes with a fixed perimeter of 1 meter.
    """
    L = 1.0  # Total length of the cut (perimeter) in meters

    print(f"--- Comparing Areas for a Fixed Perimeter of {L} meter ---\n")

    # --- Equilateral Triangle ---
    s_tri = L / 3
    area_tri = (math.sqrt(3) / 4) * s_tri**2
    print("Equilateral Triangle:")
    print(f"Equation: Area = (sqrt(3)/4) * (Perimeter/3)^2")
    print(f"Calculation: Area = (sqrt(3)/4) * ({L:.1f}/3)^2 = {area_tri:.6f} sq. meters\n")

    # --- Square ---
    s_sq = L / 4
    area_sq = s_sq**2
    print("Square:")
    print(f"Equation: Area = (Perimeter/4)^2")
    print(f"Calculation: Area = ({L:.1f}/4)^2 = {area_sq:.6f} sq. meters\n")

    # --- Regular Hexagon ---
    s_hex = L / 6
    area_hex = (3 * math.sqrt(3) / 2) * s_hex**2
    print("Regular Hexagon:")
    print(f"Equation: Area = (3*sqrt(3)/2) * (Perimeter/6)^2")
    print(f"Calculation: Area = (3*sqrt(3)/2) * ({L:.1f}/6)^2 = {area_hex:.6f} sq. meters\n")

    # --- Circle ---
    # The optimal shape according to the isoperimetric theorem.
    r_circ = L / (2 * math.pi)
    area_circ = math.pi * r_circ**2
    print("Circle (Optimal Shape):")
    print(f"Equation: Area = pi * (Perimeter / (2*pi))^2 = Perimeter^2 / (4*pi)")
    print(f"Calculation: Area = {L:.1f}^2 / (4 * {math.pi:.6f}) = {area_circ:.6f} sq. meters\n")

    print("--- Conclusion ---")
    print("The circle encloses the maximum area for a given perimeter.")
    print("The optimal cut has the symmetry of a circle.")

solve_isoperimetric_problem()