import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus of points on a cube's surface and expresses it as a percentage.
    """
    # The side length 'r' will cancel out, so we can use r=1 for calculation.
    r = 1.0
    pi = math.pi

    # --- Part 1: Locus on the two faces adjacent to P's edge ---
    # Unfolding these two faces creates a 2r x r rectangle. P is at the center of a 2r side.
    # In a local coordinate system with P at (0, r/2), the locus is part of a circle
    # x^2 + (y - r/2)^2 = r^2. The rectangle covers y in [0, r].
    # This restricts sin(theta) to [-1/2, 1/2], giving two arcs, each spanning pi/3.
    # The total angle is 2*pi/3.
    length_adjacent = r * (2 * pi / 3)

    # --- Part 2: Locus on the four "next-over" faces ---
    # This involves a path crossing an intermediate face.
    # By unfolding, we find P at a corner of a 2r x r rectangle, e.g., at (r/2, 0).
    # The locus is part of the circle (x - r/2)^2 + y^2 = r^2.
    # The arc on each of these four faces subtends an angle of pi/3.
    # The length of each arc is r * (pi/3).
    length_further = 4 * r * (pi / 3)

    # --- Part 3: Total Length and Final Calculation ---
    # Sum the lengths from all arcs.
    total_length_c = length_adjacent + length_further
    
    # The problem asks to divide this length by 2*pi*r
    divisor = 2 * pi * r
    
    # Calculate the ratio
    ratio = total_length_c / divisor
    
    # Convert the ratio to a whole number percentage
    percentage = int(round(ratio * 100))
    
    print("Let the side length of the cube be r.")
    print("The problem asks for the ratio of the length of a curve C to the value 2 * pi * r.")
    print("\nStep 1: The length of the curve on the two faces adjacent to point P is (2/3) * pi * r.")
    print(f"Step 2: The length of the curve on the four other reachable faces is (4/3) * pi * r.")
    print(f"\nStep 3: The total length of C is the sum of these parts.")
    print(f"Total Length = (2/3)*pi*r + (4/3)*pi*r = (6/3)*pi*r = 2*pi*r")
    
    print("\nFinal Equation:")
    # Here we output the numbers in the final equation using r=1 as an example.
    print(f"Length of C = {total_length_c:.4f}")
    print(f"Value to divide by = {divisor:.4f}")
    print(f"({total_length_c:.4f}) / ({divisor:.4f}) = {ratio:.4f}")

    print(f"\nExpressed as a whole number percentage, the final answer is: {percentage}%")


solve_cube_locus_problem()