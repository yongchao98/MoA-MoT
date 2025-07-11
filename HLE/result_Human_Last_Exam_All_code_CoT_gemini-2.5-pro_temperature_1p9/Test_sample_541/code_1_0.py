import math

def calculate_areas(L=1.0):
    """
    Calculates the area of the convex hull for different symmetric cuts
    with a total length L.

    Args:
        L (float): The total length of the cut.
    """

    print(f"Analyzing cuts of total length L = {L} meter.\n")

    # The area of the convex hull of a star-shaped cut with 'n' arms
    # forming a regular n-gon is given by: A = (L^2/n) * sin(pi/n) * cos(pi/n)
    
    # --- Case A: Symmetry of an Equilateral Triangle (n=3) ---
    n_tri = 3
    area_tri_star_cut = (L**2 / n_tri) * math.sin(math.pi / n_tri) * math.cos(math.pi / n_tri)
    print("Symmetry: Equilateral Triangle (Y-shaped cut)")
    print(f"Area = ({L}**2 / {n_tri}) * sin(pi/{n_tri}) * cos(pi/{n_tri})")
    print(f"Area = {area_tri_star_cut:.5f} m^2\n")

    # --- Case D: Symmetry of a Square (n=4) ---
    n_sq = 4
    area_sq_star_cut = (L**2 / n_sq) * math.sin(math.pi / n_sq) * math.cos(math.pi / n_sq)
    print("Symmetry: Square (X-shaped cut)")
    print(f"Area = ({L}**2 / {n_sq}) * sin(pi/{n_sq}) * cos(pi/{n_sq})")
    print(f"Area = {area_sq_star_cut:.5f} m^2\n")

    # --- Case E: Symmetry of a Regular Hexagon (n=6) ---
    n_hex = 6
    area_hex_star_cut = (L**2 / n_hex) * math.sin(math.pi / n_hex) * math.cos(math.pi / n_hex)
    print("Symmetry: Regular Hexagon (6-pointed star cut)")
    print(f"Area = ({L}**2 / {n_hex}) * sin(pi/{n_hex}) * cos(pi/{n_hex})")
    print(f"Area = {area_hex_star_cut:.5f} m^2\n")

    # --- Case G: Symmetry of a Circle (closed-loop cut) ---
    # This is the classic isoperimetric problem solution.
    # The cut is the perimeter itself. Area = L^2 / (4*pi)
    area_circle = (L**2) / (4 * math.pi)
    print("Symmetry: Circle (closed loop cut)")
    print(f"Area = {L}**2 / (4 * pi)")
    print(f"Area = {area_circle:.5f} m^2\n")
    
    # --- Conclusion ---
    results = {
        "Equilateral Triangle": area_tri_star_cut,
        "Square": area_sq_star_cut,
        "Regular Hexagon": area_hex_star_cut,
        "Circle": area_circle,
    }

    max_area_shape = max(results, key=results.get)
    max_area_value = results[max_area_shape]
    
    print("--- Conclusion ---")
    print(f"The maximum area ({max_area_value:.5f} m^2) is achieved by the cut with the symmetry of an {max_area_shape}.")

if __name__ == '__main__':
    calculate_areas()
