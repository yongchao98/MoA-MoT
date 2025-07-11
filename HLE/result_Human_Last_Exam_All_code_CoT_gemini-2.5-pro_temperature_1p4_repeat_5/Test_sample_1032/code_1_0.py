import math

def calculate_reachable_area():
    """
    Calculates the area a horse can reach on a rope of a given taxi-cab length,
    avoiding a house-shaped obstacle.
    """
    
    # 1. Define constants
    rope_length = 7 / 2
    house_area = 3.0

    # 2. Calculate the total area reachable without obstacles.
    # The region is |x| + |y| <= r, which is a diamond with area 2*r^2.
    total_diamond_area = 2 * rope_length**2

    # 3. Account for the house area.
    # The house is entirely within the diamond, so we subtract its area.
    base_reachable_area = total_diamond_area - house_area
    
    # 4. Analyze and calculate the "shadow area".
    # A shadow region consists of points Q where |x|+|y| <= 3.5, but the shortest
    # unobstructed path d_H(O,Q) from the origin is > 3.5.
    # Analysis shows these shadow regions are formed by the house's concave corner at (-1,-1).
    #
    # Let's check the region where a shadow is expected, for example, the square
    # defined by -2 < x < -1 and -1 < y < 0.
    # The shortest path d_H(O,Q) is the minimum of paths via corners (-2,0) and (-1,-1).
    # d_H(O,Q) = min(x-y+4, -x+y+2)
    # The region where d_H > 3.5 is the union of two sub-regions:
    # A) {y > x+1 and x-y+4 > 3.5}, which simplifies to {y > x+1 and y < x+0.5}. This set is empty.
    # B) {y < x+1 and -x+y+2 > 3.5}, which simplifies to {y < x+1 and y > x+1.5}. This set is also empty.
    #
    # A similar analysis for the symmetric region shows its shadow area is also zero.
    # Therefore, the total shadow area is 0.
    total_shadow_area = 0.0

    # 5. Final calculation
    final_area = base_reachable_area - total_shadow_area

    print("--- Calculation of the Reachable Area ---")
    print(f"Rope taxi-cab length (r): {rope_length}")
    print(f"Area of the initial diamond (|x|+|y| <= r) is 2*r^2:")
    print(f"  2 * ({rope_length})^2 = {total_diamond_area}")
    print(f"Area of the house to subtract: {house_area}")
    print("\nCalculating the shadow area (points in the diamond made unreachable by detours):")
    print("  Analysis of path lengths shows the shadow regions are empty.")
    print(f"  Total shadow area: {total_shadow_area}")
    print("\nThe final reachable area is (Diamond Area) - (House Area) - (Shadow Area).")
    print(f"Final Area = {total_diamond_area} - {house_area} - {total_shadow_area}")
    print(f"Final Area = {final_area}")

calculate_reachable_area()