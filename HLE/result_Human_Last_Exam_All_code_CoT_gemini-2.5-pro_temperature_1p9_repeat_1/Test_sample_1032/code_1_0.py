import math

def solve_horse_area():
    """
    Calculates the area a horse can reach, tied to the origin with a taxi-cab rope,
    avoiding a house.
    """
    
    # Rope length
    L = 7 / 2
    
    # Step 1: Calculate the area of the reachable diamond without any obstacles.
    # The region is |x| + |y| <= L, which is a diamond shape.
    # The area of this diamond is 2 * L^2.
    total_potential_area = 2 * L**2
    
    # Step 2: Define the area of the house obstacle.
    # The problem states the house consists of three unit squares.
    house_area = 3.0
    
    # Step 3: Justification for the calculation method.
    # The length of the shortest path from the origin (0,0) to any point (x,y)
    # around the given L-shaped house obstacle (which has a corner at the origin)
    # in the taxi-cab metric is still |x| + |y|.
    # This means a point (x,y) is reachable if and only if |x| + |y| <= L.
    # Therefore, the total reachable area is the area of the diamond minus the
    # area of the house that the horse cannot enter.
    
    # Step 4: Calculate the final area.
    final_area = total_potential_area - house_area
    
    # Step 5: Print the explanation and the final equation with numerical values.
    print(f"The horse is on a rope of taxi-cab length L = {L}.")
    print("Without any obstacles, the reachable region is a diamond shape defined by |x| + |y| <= L.")
    print(f"The area of this diamond is 2 * L^2 = 2 * {L}^2 = {total_potential_area}.")
    print(f"The horse cannot enter the house, which has an area of {house_area} unit squares.")
    print("For this specific problem, the shortest path from the origin to any point (x,y) around the house is still |x| + |y|.")
    print("Thus, the final reachable area is the diamond's area minus the house's area.")
    print(f"Final Area = {total_potential_area} - {house_area} = {final_area}")

solve_horse_area()