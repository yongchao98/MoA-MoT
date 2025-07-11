import math

def solve_triangle_grid_problem():
    """
    Calculates the maximum number of grid squares a specific triangle's perimeter can pass through.
    """
    # Step 1: Define triangle properties and optimal orientation.
    # The triangle is an isosceles right triangle with legs of length 18.
    # To maximize grid crossings, we orient the legs along slopes 2 and -1/2.
    # Let the right-angle vertex B be placed near the origin.
    epsilon = 0.001
    delta = 0.001
    B_x, B_y = epsilon, delta

    # Step 2: Calculate the coordinates of the other vertices, A and C.
    # The vector for leg BC has length 18 and slope 2.
    # The vector for leg BA has length 18 and slope -1/2.
    sqrt5 = math.sqrt(5)
    
    # Vector v_BC = (18 * 1/sqrt(5), 18 * 2/sqrt(5))
    C_x = B_x + 18 / sqrt5
    C_y = B_y + 36 / sqrt5
    
    # Vector v_BA = (18 * -2/sqrt(5), 18 * 1/sqrt(5))
    A_x = B_x - 36 / sqrt5
    A_y = B_y + 18 / sqrt5

    # Step 3: Calculate the number of squares crossed by each side.
    # Formula: N = 1 + |floor(x2) - floor(x1)| + |floor(y2) - floor(y1)|
    N_BA = 1 + abs(math.floor(A_x) - math.floor(B_x)) + abs(math.floor(A_y) - math.floor(B_y))
    N_BC = 1 + abs(math.floor(C_x) - math.floor(B_x)) + abs(math.floor(C_y) - math.floor(B_y))
    N_AC = 1 + abs(math.floor(C_x) - math.floor(A_x)) + abs(math.floor(C_y) - math.floor(A_y))

    # Step 4: Calculate the total number of unique squares.
    # The total number of squares is the size of the union of squares for each side.
    # The sets of squares overlap only at the vertex squares because the paths diverge.
    # k = N_BA + N_BC + N_AC - (overlaps at A, B, and C)
    overlap_A = 1
    overlap_B = 1
    overlap_C = 1
    total_k = N_BA + N_BC + N_AC - overlap_A - overlap_B - overlap_C

    # Print the detailed calculation.
    print("The triangle has side lengths 18, 18, and 18*sqrt(2).")
    print("To maximize the number of crossed squares, we orient its legs with slopes 2 and -1/2.")
    print("\nLet the vertices be A, B, C with the right angle at B, placed near the origin.")
    print(f"Calculated integer parts of coordinates relative to B(0,0):")
    print(f"  floor(A_x-B_x) = {math.floor(A_x - B_x)}, floor(A_y-B_y) = {math.floor(A_y - B_y)}")
    print(f"  floor(C_x-B_x) = {math.floor(C_x - B_x)}, floor(C_y-B_y) = {math.floor(C_y - B_y)}")
    
    print("\nNumber of squares crossed by each side:")
    print(f"Side BA: N_BA = 1 + |{math.floor(A_x)} - {math.floor(B_x)}| + |{math.floor(A_y)} - {math.floor(B_y)}| = {N_BA}")
    print(f"Side BC: N_BC = 1 + |{math.floor(C_x)} - {math.floor(B_x)}| + |{math.floor(C_y)} - {math.floor(B_y)}| = {N_BC}")
    print(f"Side AC: N_AC = 1 + |{math.floor(C_x)} - {math.floor(A_x)}| + |{math.floor(C_y)} - {math.floor(A_y)}| = {N_AC}")
    
    print("\nThe total number of unique squares (k) is the sum of squares for each side minus the overlaps at the vertices.")
    print(f"Final Equation: k = {N_BA} + {N_BC} + {N_AC} - 1 (at A) - 1 (at B) - 1 (at C)")
    print(f"Result: k = {total_k}")

solve_triangle_grid_problem()