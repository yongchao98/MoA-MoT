import math

def solve_cube_locus_percentage():
    """
    Solves the cube surface locus problem and returns the final answer as a percentage.
    """
    # We can assume a side length r=1.0, as the 'r' term will cancel out in the final ratio.
    r = 1.0
    
    print("### Solving the Cube Surface Locus Problem ###\n")
    print(f"The cube has a side length r = {r}.")
    print(f"We are finding the length of the locus of points on the surface at a distance r from a midpoint of an edge.\n")

    # --- Part 1: Arcs on faces containing the point P's edge ---
    # These are the two faces that share the edge P is on. By unfolding them,
    # the locus on each forms a circular arc of radius r.
    # Geometric analysis shows the angle for this arc is pi/3.
    num_arcs_part1 = 2
    arc_angle = math.pi / 3
    arc_length = r * arc_angle
    total_length_part1 = num_arcs_part1 * arc_length

    print("--- Step 1: Arcs on the two faces containing P's edge ---")
    print(f"Number of arcs: {num_arcs_part1}")
    print(f"Length of each arc = r * angle = {r} * (pi/3) = {arc_length:.4f}")
    print(f"Subtotal length: {num_arcs_part1} * {arc_length:.4f} = {total_length_part1:.4f}\n")

    # --- Part 2: Arcs on faces adjacent to the first two faces ---
    # These are formed by paths from P crossing an edge onto an adjacent face.
    # There are 4 such arcs, and by symmetry, they are identical to the first set.
    num_arcs_part2 = 4
    total_length_part2 = num_arcs_part2 * arc_length

    print("--- Step 2: Arcs on the four faces adjacent to P's edge ---")
    print(f"Number of arcs: {num_arcs_part2}")
    print(f"Length of each arc = r * angle = {r} * (pi/3) = {arc_length:.4f}")
    print(f"Subtotal length: {num_arcs_part2} * {arc_length:.4f} = {total_length_part2:.4f}\n")
    
    # --- Part 3: Total Length ---
    # The total length is the sum of all 6 arcs.
    total_length_C = total_length_part1 + total_length_part2
    
    print("--- Step 3: Total Length of Locus C ---")
    print(f"Total number of arcs = {num_arcs_part1 + num_arcs_part2}")
    print(f"Total length = 6 * (r * pi/3) = 2 * pi * r")
    print(f"Calculated Total Length: {total_length_C:.4f}\n")

    # --- Step 4: Final Ratio and Percentage ---
    # The final step is to calculate (Length of C) / (2 * pi * r).
    denominator = 2 * math.pi * r
    ratio = total_length_C / denominator
    percentage = int(round(ratio * 100))
    
    print("--- Step 4: Final Calculation ---")
    print("The final calculation is (Length of C) / (2 * pi * r)")
    # The prompt requests to output each number in the final equation.
    print(f"Final Equation: {total_length_C} / ( {2} * {math.pi} * {r} )")
    print(f"Resulting ratio: {ratio}")
    print(f"Result as a percentage: {ratio * 100}%")
    print(f"\nThe answer as a whole number percentage is: {percentage}")

solve_cube_locus_percentage()

# The final answer as a whole number.
print("\n<<<100>>>")