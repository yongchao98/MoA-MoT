import math

def find_optimal_shape():
    """
    Calculates the area of the convex hull for a star-shaped cut structure
    with N spokes and a total length of 1 meter. It demonstrates that
    the area is maximized for N=3.
    """
    L = 1.0  # Total length of the cuts

    print("--- Analysis of Convex Hull Area for a Star-Shaped Cut Structure ---")
    print(f"Total cut length L = {L} meter.")
    print("The structure consists of N spokes of length L/N radiating from a center.")
    print("The area of the convex hull is given by: Area(N) = (L^2 / (2*N)) * sin(2*pi/N)\n")

    max_area = 0.0
    optimal_n = 0

    print("N (spokes) | Area (m^2)")
    print("-----------|-----------")
    # A polygon requires at least 3 vertices
    for n in range(3, 8):
        area = (L**2 / (2 * n)) * math.sin(2 * math.pi / n)
        print(f"{n:<10} | {area:.6f}")
        if area > max_area:
            max_area = area
            optimal_n = n

    print("-" * 27)
    print(f"The maximum area is achieved when N = {optimal_n}.")
    print("This optimal cut structure has the symmetry of an equilateral triangle.")

    print("\n--- Final Equation for the Optimal Case (N=3) ---")
    N = 3
    L_val = 1.0
    # Breaking down the formula: Area = (L² / (2 * N)) * sin(2 * pi / N)
    term_L_sq = L_val**2
    term_2N = 2 * N
    term_angle_rad = (2 * math.pi) / N
    term_sin = math.sin(term_angle_rad)
    final_area = (term_L_sq / term_2N) * term_sin

    print(f"The final equation is: ({L_val}**2 / (2 * {N})) * sin(2 * pi / {N})")
    print(f"Which calculates to: ({term_L_sq} / {term_2N}) * {term_sin:.6f} = {final_area:.6f}")

find_optimal_shape()