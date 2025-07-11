import math

def calculate_locus_percentage():
    """
    Solves the cube locus problem by summing the lengths of all arcs forming the locus.
    """
    
    # Let the side length of the cube be 'r'. The final answer is independent of r's value.
    # We can perform the calculation with a placeholder value, e.g., r=1.
    r = 1.0

    # Step 1: Calculate the length of the locus C.
    # Based on geometric analysis by unfolding the cube, the locus C is found to be
    # the union of 6 distinct circular arcs on the cube's surface.
    # Each of these arcs has a calculated length of (pi * r / 3).

    # The 6 arcs are located on:
    # 1. The 'Right' face (adjacent to P).
    # 2. The 'Bottom' face (adjacent to P).
    # 3. The 'Front' face (path via the 'Right' face).
    # 4. The 'Front' face (path via the 'Bottom' face).
    # 5. The 'Back' face (path via the 'Right' face).
    # 6. The 'Back' face (path via the 'Bottom' face).
    # Loci on the 'Top' and 'Left' faces have zero length as the minimum distance is r.
    
    num_arcs = 6
    arc_length = (math.pi * r) / 3

    # Calculate the total length of the locus C.
    total_length_C = num_arcs * arc_length

    # Step 2: Define the value to divide by, which is (2 * pi * r).
    denominator = 2 * math.pi * r

    # Step 3: Compute the ratio of the locus length to the specified value.
    ratio = total_length_C / denominator

    # Step 4: Convert the ratio to a whole number percentage.
    percentage = round(ratio * 100)

    # --- Final Output ---
    # The instructions require printing each number in the final equation.
    
    print("--- Calculation Breakdown ---")
    print(f"For a cube of side length r, the locus of points is composed of {num_arcs} arcs.")
    print(f"The length of each arc is (pi * r / 3).")
    print(f"Equation for Total Length (L_C): L_C = Number of Arcs * Arc Length")
    print(f"L_C = {num_arcs} * (pi * r / 3) = {total_length_C / (math.pi * r):.2f} * pi * r")
    
    print("\nEquation for the final ratio:")
    print(f"Ratio = Total Length / (2 * pi * r)")
    print(f"{ratio:.2f} = ({total_length_C / (math.pi * r):.2f} * pi * r) / (2 * pi * r)")

    print(f"\nFinal conversion to percentage:")
    print(f"Percentage = Ratio * 100")
    print(f"Percentage = {ratio:.2f} * 100 = {percentage}%")


calculate_locus_percentage()