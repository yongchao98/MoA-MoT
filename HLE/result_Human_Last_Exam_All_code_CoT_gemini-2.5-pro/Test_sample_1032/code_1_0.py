def solve_horse_area():
    """
    Calculates the area a horse can reach, tied to the origin with a rope of a given
    taxi-cab length, with a house-shaped obstacle.
    """

    # Step 1: Define the parameters from the problem statement.
    # The rope has a taxi-cab length of 7/2.
    L = 7 / 2

    # Step 2: Calculate the total area the horse could reach without any obstacles.
    # The region is defined by |x| + |y| <= L. This is a diamond shape centered at the origin.
    # The area of this diamond is given by the formula 2 * L^2.
    total_area = 2 * L**2

    # Step 3: Define the area of the obstacle.
    # The house consists of three unit squares.
    house_area = 3

    # Step 4: Analyze the effect of the obstacle on the taxi-cab distance.
    # The horse is tied at the origin (0,0). The house is an L-shaped region in the
    # third quadrant, also touching the origin.
    # The reachable region is the set of points P(x,y) where the shortest taxi-cab
    # path from the origin to P that avoids the interior of the house has a length
    # of at most L.
    #
    # We can show that for any point P(x,y) outside the house, the shortest path
    # around the house has the same length as the standard taxi-cab distance, |x| + |y|.
    # This is because the house is made of axis-aligned squares and is located at the origin.
    # A path can always be constructed along the boundary of the house without adding
    # extra length compared to a path in an empty field.
    # For example, for a point P(x,y) with x < -1 and y < -1, a valid path is
    # (0,0) -> (-1,0) -> (-1,-1) -> (x,-1) -> (x,y).
    # The length of this path is |-1-0| + |-1-0| + |x-(-1)| + |y-(-1)| = 1 + 1 - (x+1) - (y+1) = -x - y.
    # This is the standard taxi-cab distance for a point in the third quadrant.
    #
    # Since the effective distance to any point is not increased, the reachable region
    # is simply the total area of the diamond minus the area of the house itself.

    # Step 5: Calculate the final area.
    final_area = total_area - house_area

    # Print the explanation and the final calculation.
    print(f"The taxi-cab length of the rope is L = 7/2 = {L}.")
    print("\nFirst, we calculate the total area the horse could reach without the house obstacle.")
    print(f"This area is for the region where |x| + |y| <= {L}, which forms a diamond shape.")
    print(f"The area of this diamond is 2 * L^2.")
    print(f"Total Area = 2 * ({L})^2 = {total_area}")
    print("\nNext, we consider the house, which is an obstacle of three unit squares.")
    print(f"The area of the house is {house_area} square units.")
    print("\nThe presence of the house prevents the horse from entering it. We must determine if it also blocks access to areas behind it (creates a 'shadow').")
    print("In taxi-cab geometry, because the house is made of axis-aligned squares located at the origin, a path can always be found along the house's boundary. The length of this shortest path around the house to any point (x,y) is still |x| + |y|.")
    print("Therefore, the house does not create a shadow region; it only subtracts its own area from the total reachable area.")
    print("\nThe final area is the total area minus the area of the house.")
    print(f"Final Area = Total Area - House Area")
    print(f"Final Area = {total_area} - {house_area} = {final_area}")

solve_horse_area()
<<<21.5>>>