import math

def calculate_and_print_distances():
    """
    Calculates and prints the possible normalized distances between hard spheres
    on a 2D plane for r <= 3.
    """
    print("For hard spheres on a 2D plane, the possible distances are determined by the geometry of their packing.")
    print("The most probable short-range arrangements resemble square and hexagonal lattices.")
    print("The normalized distance r = D/d, where D is the center-to-center distance and d is the sphere diameter.")
    print("\nThe set of possible distances for r <= 3 is derived below, showing the calculation for each distance:")

    # Using a dictionary to store unique distances (r_squared) and their calculation formula
    # This automatically handles duplicates from different lattice considerations.
    distances = {}

    # --- Square Lattice Calculations ---
    # r^2 = m^2 + n^2. We need m^2 + n^2 <= 9.
    for m in range(4):
        for n in range(m + 1):
            if m == 0 and n == 0:
                continue
            r_squared = m**2 + n**2
            if r_squared <= 9:
                if r_squared not in distances:
                    formula = f"sqrt({m}^2 + {n}^2)"
                    distances[r_squared] = ("square lattice", formula)

    # --- Hexagonal Lattice Calculations ---
    # r^2 = m^2 + mn + n^2. We need m^2 + mn + n^2 <= 9.
    for m in range(4):
        # We need to check negative n as well for hexagonal lattice vectors
        for n in range(-m, m + 1):
            if m == 0 and n == 0:
                continue
            r_squared = m**2 + m*n + n**2
            if r_squared <= 9:
                if r_squared not in distances:
                    # Use absolute values for cleaner formula display
                    formula = f"sqrt({abs(m)}^2 + {abs(m)}*{abs(n)} + {abs(n)}^2)" if m*n >= 0 else f"sqrt({m}^2 + {m}*({n}) + ({n})^2)"
                    distances[r_squared] = ("hexagonal lattice", formula)

    # --- Sort and Print Results ---
    sorted_r_squared = sorted(distances.keys())
    
    final_values = []
    print("\n--- Set of Distances (r <= 3) ---")
    for r_sq in sorted_r_squared:
        r_val = math.sqrt(r_sq)
        origin, eq = distances[r_sq]
        print(f"r = {r_val:.2f} (from {origin}: {eq} = sqrt({r_sq}))")
        final_values.append(f"{r_val:.2f}")
    
    # The final answer format as requested by the prompt
    print("\n<<<" + ", ".join(final_values) + ">>>")

calculate_and_print_distances()