import math

def get_area_formula():
    """
    This function derives and prints the formula for the area of the triangle T(t).
    """

    # Constants from the problem
    R = 10  # Radius of the circle
    S = 10  # Side length of the hexagon (S=R)
    v = 1   # Speed of the vertices of the triangle

    # --- Step 1: Explain the reasoning ---
    print("Step 1: The rotation of the hexagon does not change the area of the inscribed triangle.")
    print("We can analyze the area in the hexagon's own rotating frame.\n")

    print("Step 2: The triangle T(t) remains equilateral due to symmetry.")
    print("Its area A(t) depends only on its side length L(t): A(t) = (sqrt(3)/4) * L(t)^2.\n")

    print("Step 3: Let s(t) be the distance of a triangle vertex from the midpoint of the hexagon side it moves on.")
    print("Through geometric analysis, the square of the triangle's side length is found to be:")
    num_1 = 225
    num_2 = 3
    print(f"L(t)^2 = {num_1} + {num_2} * s(t)^2\n")

    print("Step 4: The function s(t) describes the oscillatory motion of a vertex on a side of length 10.")
    print("The motion is periodic with a period of 20 seconds and an amplitude of 5 units.")
    print("s(t) is defined piece-wise for t' = t % 20:")
    print("  s(t) = t'          if 0 <= t' <= 5")
    print("  s(t) = 10 - t'     if 5 < t' <= 15")
    print("  s(t) = t' - 20     if 15 < t' <= 20\n")

    # --- Step 5: Print the final combined formula ---
    print("Step 5: Combining these results gives the final equation for the area A(t).")
    print("The final equation is:")
    
    final_num_sqrt = 3
    final_num_div = 4
    final_num_add = 225
    final_num_mult = 3

    # Outputting each number in the final equation as requested.
    print(f"A(t) = (sqrt({final_num_sqrt}) / {final_num_div}) * ({final_num_add} + {final_num_mult} * s(t)^2)")
    print("where s(t) is the periodic function defined in Step 4.")

# Execute the function to print the solution
get_area_formula()