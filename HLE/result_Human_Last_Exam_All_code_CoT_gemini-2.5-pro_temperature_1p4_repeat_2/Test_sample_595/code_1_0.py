import math

def calculate_k():
    """
    Calculates the largest number k of grid squares the triangle's perimeter can pass through.
    """
    
    # Side lengths
    side_a = 18
    side_b = 18
    side_c = 18 * math.sqrt(2)
    
    # The triangle is a right-angled isosceles triangle.
    # To maximize grid crossings, we align the legs of the right angle with the diagonals of the grid.
    # The vectors for the two legs of length 18 are v1 and v2.
    # To maximize Manhattan distance |x|+|y| for a given length L, we need x=y=L/sqrt(2).
    # For length 18, the components are 18/sqrt(2) = 9*sqrt(2).
    
    # The vectors representing the sides from the right-angled vertex A are:
    # v_AB = (9*sqrt(2), 9*sqrt(2))
    # v_AC = (9*sqrt(2), -9*sqrt(2))
    
    dx_ab = 9 * math.sqrt(2)
    dy_ab = 9 * math.sqrt(2)
    
    dx_ac = 9 * math.sqrt(2)
    dy_ac = -9 * math.sqrt(2)
    
    # Let the triangle vertices be A, B, C.
    # Let the floor-projected vertices be FA, FB, FC.
    # We place FA at the origin of the integer grid.
    FA_x, FA_y = 0, 0
    
    # The position of FB and FC depends on the fractional part of vertex A's coordinates.
    # Let A_real = (f_x, f_y) where f_x, f_y are in [0, 1).
    # FB = (floor(f_x + dx_ab), floor(f_y + dy_ab))
    # FC = (floor(f_x + dx_ac), floor(f_y + dy_ac))
    # We want to choose f_x and f_y to maximize the Manhattan perimeter of the triangle (FA, FB, FC).
    
    # To maximize floor(f + val), we choose f such that f + frac(val) >= 1.
    # frac(9*sqrt(2)) = frac(12.7279...) = 0.7279...
    # We want f_x > 1 - 0.7279 = 0.2721
    # We want f_y > 1 - 0.7279 = 0.2721
    # Let's choose f_x = 0.3 and f_y = 0.8 to satisfy this and the non-lattice point conditions.
    
    f_x = 0.3
    f_y = 0.8
    
    # Coordinates of the floor-projected vertices relative to FA.
    FB_x = math.floor(f_x + dx_ab)
    FB_y = math.floor(f_y + dy_ab)
    
    FC_x = math.floor(f_x + dx_ac)
    FC_y = math.floor(f_y + dy_ac)

    # Now, calculate the Manhattan perimeter of the integer triangle (FA, FB, FC)
    
    # Distance from FA(0,0) to FB
    dist_ab = abs(FB_x - FA_x) + abs(FB_y - FA_y)
    
    # Distance from FB to FC
    dist_bc = abs(FC_x - FB_x) + abs(FC_y - FB_y)
    
    # Distance from FC to FA(0,0)
    dist_ca = abs(FA_x - FC_x) + abs(FA_y - FC_y)
    
    k = dist_ab + dist_bc + dist_ca
    
    print(f"The side lengths of the triangle are {side_a}, {side_b}, and {side_c:.2f}.")
    print("To maximize grid crossings, the legs of the right triangle are oriented diagonally.")
    print("The displacement vectors for the legs are approximately ({dx_ab:.2f}, {dy_ab:.2f}) and ({dx_ac:.2f}, {dy_ac:.2f}).")
    print(f"By choosing an optimal starting position for vertex A, the integer grid vertices are:")
    print(f"F(A) = ({FA_x}, {FA_y})")
    print(f"F(B) = ({FB_x}, {FB_y})")
    print(f"F(C) = ({FC_x}, {FC_y})")
    print("\nThe number of squares crossed is the Manhattan perimeter of triangle F(A)F(B)F(C).")
    print(f"Number of squares for side AB = |{FB_x} - {FA_x}| + |{FB_y} - {FA_y}| = {dist_ab}")
    print(f"Number of squares for side BC = |{FC_x} - {FB_x}| + |{FC_y} - {FB_y}| = {dist_bc}")
    print(f"Number of squares for side CA = |{FA_x} - {FC_x}| + |{FA_y} - {FC_y}| = {dist_ca}")
    print(f"\nTotal squares k = {dist_ab} + {dist_bc} + {dist_ca} = {k}")
    
    # Return the final value of k.
    return k

final_k = calculate_k()
print(f"\nThe largest number k is: {final_k}")
<<<76>>>